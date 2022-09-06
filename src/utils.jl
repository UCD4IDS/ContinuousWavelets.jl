"""
    morsefreq(c::CWT{W,T,Morse,N}) where {W,T,N}

Measures of frequency for generalized Morse wavelet. [with F. Rekibi]
The output returns the modal or peak.
For `be=0`, the wavelet becomes an analytic lowpass filter.
Reference: Lilly and Olhede (2009).  Higher-order properties of analytic wavelets.
IEEE Trans. Sig. Proc., 57 (1), 146--160.
"""
function morsefreq(c::CWT{W,T,Morse,N}) where {W,T,N}
    ga = c.waveType.ga
    be = c.waveType.be

    fm = @. exp((log(be) - log(ga)) / ga)

    if sum(be .== 0) != 0 && isempty(fm)
        fm = (log(2))^(1 / ga)
    elseif sum(be .== 0) != 0 && !isempty(fm)
        fm[be.==0] = (log(2))^(1 / ga[be.==0])
    end

    return fm
end


"""
    getNWavelets(n1,c) -> nOctaves, totalWavelets, sRanges, sWidth

Utility for understanding the spacing of the wavelets. `sRanges` is a list of
the s values used in each octave. `sWidth`` is a list of the corresponding
variance adjustments.
"""
function getNWavelets(n1, c::CWT)
    nOctaves = getNOctaves(n1, c)
    n, nSpace = setn(n1, c)
    isAve = !(typeof(c.averagingType) <: NoAve)
    if nOctaves ≤ c.averagingLength + getMinScaling(c) # there isn't enough space for anything but averaging.
        totalWavelets = isAve
        sRanges = [1.0]
        sWidth = [1.0]
        return nOctaves, totalWavelets, sRanges, sWidth
    end
    sRange = 2 .^ (polySpacing(nOctaves, c))
    totalWavelets = round(Int, length(sRange) + isAve)
    sWidth = varianceAdjust(c, totalWavelets, nOctaves)
    return nOctaves, totalWavelets, sRange, sWidth
end


"""
Different wavelet familes need to end at a different number of octaves because they have different tail behavior.
"""
getNOctaves(n1, c::CWT{W,T,M,N}) where {W,T,N,M} = log2(n1 >> 1 + 1) + c.extraOctaves
# choose the number of octaves so the last mean, which is at s*σ[1]
# is 3 standard devations away from the end
getNOctaves(n1, c::CWT{W,T,Morlet,N}) where {W,T,N} = log2((n1 >> 1 + 1) / (c.σ[1] + 3)) + c.extraOctaves
getNOctaves(n1, c::CWT{W,T,<:Paul,N}) where {W,T,N} = log2((n1 >> 1 + 1) / (2c.α + 5)) + c.extraOctaves
# choose the number of octaves so the last mean is 5 standard deviations from the end
function getNOctaves(n1, c::CWT{W,T,<:Dog,N}) where {W,T,N}
    μ = getMean(c)
    σ = getStd(c)
    log2(n1 >> 1 / (μ + 5σ)) + c.extraOctaves
end
# choose the number of octaves so the smallest support is twice the qmf
getNOctaves(n1, c::CWT{W,T,<:ContOrtho,N}) where {W,T,N} = log2(n1) - 2 - log2(length(qmf(c.waveType))) + c.extraOctaves
getNOctaves(n1, c::CWT{W,T,Morse,N}) where {W,T,N} = log2((n1 >> 1 + 1) / (morsefreq(c) + 1)) + c.extraOctaves
# getNOctaves(n1,c::CWT{W,T, Morse, N}) where {W, T, N} = 4 + c.extraOctaves

"""
As with the last octave, different wavelet families have different space decay rates, and in the case of symmetric or zero padding we don't want wavelets that bleed across the boundary from the center.
"""
getMinScaling(c::CWT{W,T,M,N}) where {W,T,N,M} = 0 # by default all scales are allowed (all of the orthogonal transforms)
getMinScaling(c::CWT{W,T,<:Morlet,N}) where {W,T,N} = 1 / (c.β)^0.8 # morlet is slightly too large at the boundaries by default
getMinScaling(c::CWT{W,T,<:Paul,N}) where {W,T,N} = 2 / (2c.α + 1) # Paul presents some difficulties, as the decay changes quickly (like 1/t^(α+1))
getMinScaling(c::CWT{W,T,<:Dog,N}) where {W,T,N} = 2 # like morlet, the decay for Dog is exponential and consistent across derivatives
# getMinScaling(c::CWT{W,T,<:Morse,N}) where {W,T,N,M} = log2(2 * morsefreq(c))
getMinScaling(c::CWT{W,T,<:Morse,N}) where {W,T,N,M} = 5 * morsefreq(c)

function varianceAdjust(this::CWT{W,T,M,N}, totalWavelets, nOct) where {W,T,N,M}
    # increases the width of the wavelets by σ[i] = (1+a(total-i)ᴾ)σₛ
    # which is a polynomial of order p going from 1 at the highest frequency to
    # sqrt(p) at the lowest
    β = this.β
    x = exp(nOct) / (100 + exp(nOct)) # the fewer octaves there are, the smaller the width adjustment we need
    a = (β^x - 1) / (totalWavelets - 1)^β
    sWidth = 1 .+ a .* (totalWavelets .- (1:totalWavelets)) .^ β
    if any(isnan.(sWidth))
        return [1]
    end
    return sWidth
end


function polySpacing(nOct, c)
    a = getMinScaling(c) + c.averagingLength
    β = c.β
    Q = c.Q
    if nOct ≤ a
        # averagingLength and the min scaling are too large for anything to be done
        return [1.0]
    elseif 0 < nOct - a ≤ 1
        # there's only one octave, just return the linear rate of Q
        return range(a, nOct, length=1 + round(Int, Q))
    end
    # x is the index, y is the scale
    # y= a + b*x^(1/β), solve for b and the final point x₁ so that
    # y(x₁) = nOct
    # dy/dx(x₁) = 1/Q
    lastWavelet = c.Q * (nOct - a)
    b = 1 / c.Q * lastWavelet^((β - 1) / β)
    # step size so that there are actually Q wavelets in the last octave
    startOfLastOctave = ((nOct - a - 1) / b)^β
    stepSize = (lastWavelet - startOfLastOctave) / Q
    roundedStepSize = lastWavelet / round(Int64, lastWavelet / stepSize)
    samplePoints = range(0, stop=lastWavelet,
        step=roundedStepSize)
    return a .+ b .* (samplePoints) .^ (1 / β)
end

# a utility to just get the start, stop, and step size used in polySpacing. Only used for explanatory purposes
function genSamplePoints(nOct, c)
    a = getMinScaling(c) + c.averagingLength
    β = c.β
    Q = c.Q
    if nOct ≤ a
        # averagingLength and the min scaling are too large for anything to be done
        return [1.0]
    elseif 0 < nOct - a ≤ 1
        # there's only one octave, just return the linear rate of Q
        return range(a, nOct, length=1 + round(Int, Q))
    end
    # x is the index, y is the scale
    # y= a + b*x^(1/β), solve for b and the final point x₁ so that
    # y(x₁) = nOct
    # dy/dx(x₁) = 1/Q
    lastWavelet = c.Q * (nOct - a)
    b = 1 / c.Q * lastWavelet^((β - 1) / β)
    # step size so that there are actually Q wavelets in the last octave
    startOfLastOctave = ((nOct - a - 1) / b)^β
    stepSize = (lastWavelet - startOfLastOctave) / Q
    roundedStepSize = lastWavelet / round(Int64, lastWavelet / stepSize)
    0, lastWavelet, roundedStepSize
end

"Adjust the length of the storage based on the boundary conditions."
function setn(n1, c)
    if boundaryType(c) <: ZPBoundary
        base2 = ceil(Int, log2(n1))   # power of 2 nearest to n1
        nSpace = 2^(base2)
        n = nSpace >> 1 + 1
    elseif boundaryType(c) <: SymBoundary
        # n1+1 rather than just n1 because this is going to be used in an rfft
        # for real data
        n = n1 + 1
        nSpace = 2 * n1
    elseif boundaryType(c) <: PerBoundary
        n = n1 >> 1 + 1
        nSpace = n1
    else
        error("No such boundary as $(boundaryType(c))")
    end
    return n, nSpace
end


# computing averaging function utils

# this makes sure the averaging function is real and positive
adjust(c::CWT) = 1
adjust(c::CWT{W,T,<:Dog}) where {W,T} = 1 / (-im)^(c.α)

function locationShift(c::CWT{W,T,<:Morlet,N}, s, ω, sWidth) where {W,T,N}
    s0 = 3 * c.σ[1] * s * sWidth
    ω_shift = ω .+ c.σ[1] * s0
    return (s0, ω_shift)
end

function locationShift(c::CWT{W,T,<:Dog,N}, s, ω, sWidth) where {W,T,N}
    μNoS = getMean(c)
    μLast = μNoS * (2s)
    s0 = μLast / 2 / getStd(c)
    μ = μNoS * s0
    ω_shift = ω .+ μ
    return (s0, ω_shift)
end

function locationShift(c::CWT{W,T,<:Morse,N}, s, ω, sWidth) where {W,T,N}
    # s0 = c.waveType.cf * s * sWidth
    s0 = morsefreq(c) * s * sWidth
    # ω_shift = ω .+ c.waveType.cf * s0
    ω_shift = ω .+ s0
    return (s0, ω_shift)
end


"Get the mean of the mother wavelet where s=1."
function getMean(c::CWT{W,T,<:Dog}, s=1) where {W,T}
    Gm1 = gamma((c.α + 1) / 2)
    Gm2 = gamma((c.α + 2) / 2)
    Gm2 / Gm1 * sqrt(2) * s^2
end
getMean(c::CWT{W,T,<:Paul}, s=1) where {W,T} = (c.α + 1) * s
function getMean(c::CWT{W,T,<:Morlet}, s=1) where {W,T}
    return s * c.σ[1]
end
function getMean(c::CWT{W,T,<:Morse}, s=1) where {W,T}
    #return s*c.waveType.cf
    return s * morsefreq(c)
end

"""
    getStd(c::CWT{W, T, <:Dog}, s=1) where {W,T}

Get the standard deviation of the mother wavelet.
"""
getStd(c::CWT{W,T,<:Dog}, s=1) where {W,T} = sqrt(c.α + 1 - getMean(c)^2) * s^(3 / 2)
getStd(c::CWT{W,T,<:Paul}, s=1) where {W,T} = sqrt((c.α + 2) * (c.α + 1)) * s
getStd(c::CWT{W,T,<:Morlet}, s=1) where {W,T} = (s * c.β^0.8) / sqrt(2)
getStd(c::CWT{W,T,<:Morse}, s=1) where {W,T} = (s * c.β^0.8) / sqrt(2)

function locationShift(c::CWT{W,T,<:Paul,N}, s, ω, sWidth) where {W,T,N}
    s0 = s * sqrt(c.α + 1)
    ω_shift = ω .+ (c.α .+ 1) * s0
    return (s0, ω_shift)
end


function getUpperBound(c::CWT{W,T,<:Morlet,N}, s) where {W,T,N}
    return c.σ[1] * s
end

function getUpperBound(c::CWT{W,T,<:Dog,N}, s) where {W,T,N}
    return sqrt(2) * s * gamma((c.α + 2) / 2) / gamma((c.α + 1) / 2)
end

function getUpperBound(c::CWT{W,T,<:Paul,N}, s) where {W,T,N}
    return (c.α + 1) * s
end

"""
    getMeanFreq(Ŵ, fsample=2000) -> arrayOfFreqs

Calculate each of the mean frequencies of a collection of analytic or real wavelets `Ŵ`.
The default sampling rate `fsample=2kHz`, so the maximum frequency is 1kHz.
"""
function getMeanFreq(Ŵ::Array, fsample=2000)
    eachNorm = [norm(w, 1) for w in eachcol(Ŵ)]
    freqs = range(0, fsample / 2, length=size(Ŵ, 1))
    return map(ŵ -> sum(abs.(ŵ) .* freqs), eachcol(Ŵ)) ./ eachNorm
end

function getMeanFreq(n1, cw::CWT, fsample=2000)
    Ŵ, ω = computeWavelets(n1, cw)
    getMeanFreq(Ŵ, fsample)
end


# orthogonal wavelet utils
function getContWaveFromOrtho(c, N)
    depth, ψLen = calcDepth(c, N)
    φ, ψ, t = makewavelet(c.waveType.o, depth)
    return φ, -ψ, ψLen
end

"""
    calcDepth(w,N) -> nIters, sigLength

Given a `CWT` with an orthogonal wavelet type, calculate the depth `nIters`
necessary to get a resulting wavelet longer than `N`. `sigLength` is the
resulting length.
"""
function calcDepth(w, N)
    origLen = length(qmf(w.waveType))
    netLen = origLen
    nIters = 0
    if N < netLen
        return nIters
    end
    while N > netLen - 2^nIters + 1
        netLen = nextLength(netLen, origLen)
        nIters += 1
    end
    return (nIters, netLen - 2^nIters + 1)
end
# the actual length calculator, as this is recursive
nextLength(curr, orig) = curr * 2 + orig - 1 # upsample then convolve
# pad a vector `v` to have length `N`. Or chop off from both sides if its too long
function padTo(v, N)
    nDiff = N - length(v)
    if nDiff > 0 # there should be more entries than there are
        cat(v, zeros(eltype(v), max(N - length(v), 0)), dims=1)
    elseif nDiff == 0
        return v
    else # there are too many, so chop off half the difference from each side
        inds = ceil(Int, 1 - nDiff / 2):length(v)+ceil(Int, nDiff / 2)
        return v[inds]
    end
end


# create interpolater for the orthogonal cases
genInterp(ψ) = interpolate(ψ, BSpline(Quadratic(Reflect(OnGrid()))))




"""
    computeDualWeights(Ŵ) -> β

Compute the weight given to each wavelet so that in the Fourier domain, the sum across
wavelets is as close to 1 at every frequency below the peak of the last wavelet
(that is `β' .*abs.(Ŵ[1:lastFreq,:]) ≈ ones(lastFreq)` in an ℓ^2 sense).
"""
function computeDualWeights(Ŵ, wav)
    @views lastReasonableFreq = computeLastFreq(Ŵ[:, end], wav)
    Wdag = pinv(Ŵ[1:lastReasonableFreq, :])
    β = conj.(Wdag * ones(size(Wdag, 2)))'
    return β
end
# function computeDualWeights(Ŵ, wav::CWT{W,T,<:ContOrtho,N}, n,naive=true) where {W,T,N}
#     β = computeNaiveDualWeights(Ŵ, wav, n)
#     return β
# end
computeLastFreq(ŵ, wav) = length(ŵ) - 3 # we have constructed the wavelets so the last frequency is <1%, so trying to get it to 1 will fail
computeLastFreq(ŵ, wav::CWT{W,T,<:Morlet}) where {W,T} = findlast(abs.(ŵ) .> 1 / 4 * maximum(abs.(ŵ))) # the gaussian decay makes this way too large by the end

"""
    computeNaiveDualWeights(Ŵ, wav, n1)

Compute the dual weights using the scaling amounts and the normalization power p. Primarily
used when the least squares version is poorly constructed. When the least squares version
doesn't perform well, this is also likely to have poor reconstruction, but it won't give
extremely negative or oscillatory weights like the least squares version.
"""
function computeNaiveDualWeights(Ŵ, wav, n1)
    _, _, sRange, _ = ContinuousWavelets.getNWavelets(n1, wav)
    isAve = !(wav.averagingType isa NoAve)
    p = wav.p
    # every wavelet is multiplied by 1/s^(1/p), so we need to multiply by that to normalize in frequency. The extra +1 is to convert from uniform in scale to uniform in frequency
    sToTheP = sRange' .^ (1 / p + 1)
    if isAve && size(Ŵ, 2) > 1
        # averaging needs to have maximum value of 1
        aveMultiplier = 1 / norm(Ŵ[:, 1], Inf)
        dual = sum(Ŵ[:, 2:end] .* sToTheP, dims=2)
        renormalize = norm(dual, Inf)
        dual ./= renormalize
        #dual += Ŵ[:,1] .* sToTheP[1]
        β = cat(aveMultiplier, 2 .* sToTheP ./ renormalize..., dims=2)
    elseif size(Ŵ, 2) > 1
        dual = sum(Ŵ .* sToTheP, dims=2)
        β = sToTheP ./ norm(dual, Inf)
        β[2:end] .* 2
    else
        # there's only the averaging
        β = [1 / norm(Ŵ[:, 1], Inf)]
    end
    return β
end

"""
    getDualCoverage(n,cWav, invType) -> dualCover, dualNorm

Get the sum of the weights and its deviation from 1.
"""
function getDualCoverage(n, cWav, invType)
    Ŵ = computeWavelets(n, cWav)[1]
    if invType isa DualFrames
        dualCover = sum(conj.(Ŵ) .* explicitConstruction(Ŵ), dims=2)
        return dualCover, norm(dualCover .- 1)
    elseif invType isa PenroseDelta
        β = computeDualWeights(Ŵ, cWav)
    elseif invType isa NaiveDelta
        β = computeNaiveDualWeights(Ŵ, cWav, n)
    end
    dualCover = sum(conj.(β) .* Ŵ, dims=2)
    dualCover, norm(dualCover .- 1)
end

function dualDeviance(n, cWav, naive)
    getDualCoverage(n, cWav, naive)
end

"""
    explicitConstruction(Ŵ)

Explicitly construct the canonical dual frame elements for a translation invariant frame.
This is likely to end poorly due to the low representation at high frequencies.
"""
function explicitConstruction(Ŵ)
    Ŵdual = conj.(Ŵ ./ [norm(ŵ)^2 for ŵ in eachslice(Ŵ, dims=1)])
end



"""
    caveats(n1, c::CWT{B,CT,W}; coiTolerance = exp(-2), fsample = 2000) -> sRange, meanFreqs, coi

Given a length `n1` and a CWT struct `c`, returns the scales `sRange` used, the mean frequencies of the wavelets `meanFreqs`, and the cone of influence `coi` for each wavelet.
Returns the period, the scales, and the cone of influence for the given wavelet transform.
If you have sampling information, you will need to scale the vector scale appropriately by
1/δt, and the actual transform by δt^(1/p).
"""
function caveats(n1, c::CWT; coiTolerance=exp(-2), fsample=2000)
    nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1, c)
    # padding determines the actual number of elements
    n, nSpace = setn(n1, c)
    # indicates whether we should keep a spot for the father wavelet
    isAve = !(typeof(c.averagingType) <: NoAve)

    # Fourier equivalent frequencies
    Ŵ, ω = computeWavelets(n1, c)
    freqs = getMeanFreq(Ŵ, fsample)
    coi = directCoiComputation(n1, c; coiTolerance=coiTolerance)
    return sRange, freqs, coi
end

caveats(Y::AbstractArray{T}, w::ContWaveClass) where {T<:Number,S<:Real} = caveats(size(Y, 1), CWT(w), J1=J1)

"""
    directCoiComputation(n1, c::CWT; coiTolerance = exp(-2)) -> coi

A straightforward computation of the cone of influence directly from the
constructed wavelets. `coi` is a binary matrix of dimensions `(signalLength)×(nscales+1)`
that indicates when the autocorrelation of the wavlet is greater than
`coiTolerance` of the maximum value.
"""
function directCoiComputation(n1, c::CWT; coiTolerance=exp(-2))
    ψ, ω = ContinuousWavelets.computeWavelets(n1, c, space=true)
    autoCorr = cat([sum(ψ .* conj.(circshift(ψ, (ii, 0))), dims=1) for ii = 1:size(ψ, 1)]..., dims=1) # auto correlation in time domain
    normedAutoCorr = autoCorr ./ maximum(abs.(autoCorr), dims=1) # normalize by the max norm
    return abs.(normedAutoCorr) .>= coiTolerance # the cone is where the autocorrelation is larger than the tolerance
end

"""
    coi(n, s, wave::CWT)

for a signal of length `n` and wavelet transform of type wave, return a binary vector indicating the cone of influence from the edges so that the auto-correlation of the wavelet has decreased by a factor of ``\\textrm{e}^2``.
"""
function coi(n, s, ::Morlet)
    sqrt(2) * s
end

@doc raw"""
    crossSpectrum(X, Y, c::CWT)

Compute the cross spectrum of signals `X` and `Y`, which is defined as ``S(\\conj{C(X)})``
Note that unlike `cwt`, `crossSpectrum` only works for collections `X` and `Y` of vectors of shape (signalLength)×(nSignals), and outputs a collection of cross spectrums that has shape (signalLength)×(nscales+1)×(nSignalsX)×(nSignalsY).

# Examples
```@meta
DocTestFilters = [r"\@ ContinuousWavelets .*", r"[0-9]\.[0-9]{5}e-[0-9][5-9]"]
```
```jldoctest
julia> using ContinuousWavelets, Random

julia> rng = MersenneTwister(23425); Y = randn(rng, 2053, 4);

julia> X = Y .+ 3;

julia> c = wavelet(morl, β = 2)

julia> Xspec = crossSpectrum(X, Y, c); size(Xspec)
(2053, 29, 4, 4)

julia> Xspec[:,:,1,1]
2053×29 Matrix{ComplexF64}:
 -4.14517e-5+2.19692e-20im  …  1.19877e-5-7.07215e-15im
 -4.14157e-5+2.23209e-21im     1.19896e-5-7.06562e-15im
            ⋮               ⋱
 0.000119144+4.38332e-21im     1.70054e-5+1.85809e-15im
 0.000119178+1.3884e-20im      1.69993e-5+1.8598e-15im

julia> Xspec[:,:,1,2]
2053×29 Matrix{ComplexF64}:
  5.42995e-5-1.94343e-20im  …    2.649e-6-1.22869e-6im
   5.4303e-5-1.52994e-20im      2.6479e-6-1.23329e-6im
            ⋮               ⋱
 -3.17457e-5-1.12611e-20im     4.71683e-6+3.72814e-6im
 -3.17719e-5+1.44436e-20im     4.71417e-6+3.7279e-6im

```
"""
function crossSpectrum(X, Y, c::CWT{B,S,W,N,isAn}) where {B,S,W,N,isAn}
    aveCXY, _ = sharedCrossSpectrum(X, Y, c)
    return aveCXY
end

@doc raw"""
    waveletCoherence(X, Y, c)

Compute the wavelet coherence between `X` and `Y`. This is given by the power of the cross spectrum, normalized by smoothed powers of both `X` and `Y`. Explicitly, that is ``\frac{|S(C^*(X)C(Y))|^2}{S(|C(X)|^2)S(|C(Y)|^2)}``.

# Examples
```@meta
DocTestFilters = r"\@ ContinuousWavelets .*"
```
```jldoctest ex
julia> using ContinuousWavelets, Random

julia> rng = MersenneTwister(23425); Y = randn(rng, 2053, 4);

julia> X = Y .+ 3;

julia> c = wavelet(morl,β=2)

julia> wCo = waveletCoherence(X, Y, c);


julia> wCo[:,:,1,1]
2053×29 Matrix{Float64}:
 0.688432  1.0  1.0  1.0  1.0  …  1.0  1.0  1.0  1.0  1.0
 0.687333  1.0  1.0  1.0  1.0     1.0  1.0  1.0  1.0  1.0
 ⋮                             ⋱       ⋮
 0.952408  1.0  1.0  1.0  1.0     1.0  1.0  1.0  1.0  1.0
 0.952626  1.0  1.0  1.0  1.0     1.0  1.0  1.0  1.0  1.0

julia> wCo[:,:,1,2]
2053×29 Matrix{Float64}:
 0.995092  0.984109  0.952974  …  0.128329  0.0335945
 0.995082  0.984082  0.952852     0.128308  0.0336075
 ⋮                             ⋱
 0.652498  0.994972  0.992932     0.19849   0.0902301
 0.653111  0.994978  0.992951     0.198585  0.0902469

```
"""
function waveletCoherence(X, Y, c)
    (aveCXY, aveWavelet, CX, CY) = sharedCrossSpectrum(X, Y, c)
    aveCX = abs.(dropdims(ContinuousWavelets.cwt(abs.(CX) .^ 2, c, aveWavelet), dims=2))
    aveCY = abs.(dropdims(ContinuousWavelets.cwt(abs.(CY) .^ 2, c, aveWavelet), dims=2))
    return abs.(aveCXY) .^ 2 ./ (aveCX .* aveCY)
end

"""
`waveletCoherence` depends on results precomputed by `crossSpectrum` in addition to the actual output of `crossSpectrum`.
"""
function sharedCrossSpectrum(X, Y, c)
    # ensure that the wavelet transform uses some sort of averaging
    isAve = !(typeof(c.averagingType) <: NoAve)
    if !isAve
        cAve = CWT(c.waveType, c.Q, ContinuousWavelets.boundaryType(c)(), ContinuousWavelets.Father(), c.averagingLength, c.frameBound, c.p, c.β)
    else
        cAve = c
    end

    # make sure both X and Y have the right shape
    if ndims(X) < 2
        X = reshape(X, (size(X)..., 1))
    elseif ndims(X) > 2
        X = reshape(X, (size(X, 1), :))
    end
    if ndims(Y) < 2
        Y = reshape(Y, (size(Y)..., 1))
    elseif ndims(Y) > 2
        Y = reshape(Y, (size(Y, 1), :))
    end
    daughters, _ = ContinuousWavelets.computeWavelets(size(Y, 1), cAve)
    aveWavelet = daughters[:, 1:1] ./ norm(daughters[:, 1:1], c.p) # we need an averaging function that returns the ave value, rather than the one that preserves the frame bound
    CX = ContinuousWavelets.cwt(X, cAve, daughters)
    CY = ContinuousWavelets.cwt(Y, cAve, daughters)
    # reshape so we can do an outer product
    CX = reshape(CX, (size(CX)[1:3]..., 1))
    CY = reshape(CY, (size(CY)[1:2]..., 1, size(CY, 3)))
    CXY = conj.(CX) .* CY # conjugate of the first times the second
    aveCXY = dropdims(ContinuousWavelets.cwt(CXY, cAve, aveWavelet), dims=2) # get the average of the multiplication
    return (aveCXY, aveWavelet, CX, CY)
end
