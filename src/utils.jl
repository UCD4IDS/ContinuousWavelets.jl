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
function getMeanFreq(Ŵ, fsample=2000)
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
```jldoctest
julia> using ContinuousWavelets, Random

julia> rng = MersenneTwister(23425); Y = randn(rng, 2053, 4);

julia> X = Y .+ 3;

julia> c = wavelet(morl)

julia> Xspec = crossSpectrum(X, Y, c); size(Xspec)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.038173051785201154
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6
(2053, 18, 4, 4)

julia> Xspec[:,:,1,1]
2053×18 Matrix{ComplexF64}:
 -5.24953e-5+5.40781e-20im  …  8.04088e-6-2.27563e-14im
 -5.25855e-5+1.07093e-20im     8.04716e-6-2.25757e-14im
  -5.2766e-5-3.30814e-20im     8.06082e-6-2.2219e-14im
 -5.30372e-5-1.06163e-20im     8.08398e-6-2.16944e-14im
 -5.33999e-5+6.16927e-21im     8.11969e-6-2.10141e-14im
 -5.38549e-5+5.36968e-20im  …  8.17177e-6-2.01937e-14im
 -5.44031e-5-3.43429e-20im     8.24462e-6-1.92514e-14im
 -5.50459e-5+3.346e-20im       8.34296e-6-1.82075e-14im
 -5.57845e-5-2.23369e-21im     8.47159e-6-1.70835e-14im
 -5.66204e-5-2.33671e-21im     8.63515e-6-1.59018e-14im
            ⋮               ⋱
  5.75205e-5+7.53592e-21im     1.17436e-5+4.49371e-15im
  4.89493e-5+1.03129e-20im  …  1.11615e-5+4.78868e-15im
  4.12281e-5+2.25653e-20im     1.06369e-5+5.06263e-15im
  3.44424e-5+1.72952e-20im     1.01758e-5+5.30987e-15im
  2.86679e-5+3.46252e-20im     9.78334e-6+5.52514e-15im
  2.39688e-5+1.17558e-20im     9.46397e-6+5.70358e-15im
  2.03977e-5+1.56031e-20im  …  9.22126e-6+5.8412e-15im
  1.79944e-5+7.12807e-20im     9.05792e-6+5.93478e-15im
   1.6786e-5+2.05713e-20im     8.97579e-6+5.98214e-15im

julia> Xspec[:,:,1,2]
2053×18 Matrix{ComplexF64}:
   0.00010865+1.01396e-20im  …  4.51034e-6+5.79035e-6im
  0.000108391+4.88923e-20im     4.50665e-6+5.73936e-6im
  0.000107875-3.90038e-20im     4.49934e-6+5.63796e-6im
  0.000107105+2.30429e-20im     4.48854e-6+5.48733e-6im
  0.000106087+9.41468e-20im     4.47439e-6+5.28921e-6im
  0.000104828-6.83516e-20im  …  4.45701e-6+5.04588e-6im
  0.000103335-1.11258e-20im     4.43643e-6+4.76011e-6im
  0.000101619+1.83624e-20im     4.41256e-6+4.43514e-6im
   9.96888e-5-9.42864e-22im      4.3851e-6+4.07463e-6im
   9.75555e-5+2.22261e-21im     4.35352e-6+3.68257e-6im
             ⋮               ⋱
 -0.000158918+2.76457e-20im     5.39643e-6+2.66862e-6im
 -0.000162695+4.25767e-20im  …  5.26804e-6+2.31855e-6im
 -0.000165975-1.70655e-20im     5.13897e-6+2.00929e-6im
 -0.000168767+2.54067e-20im     5.01545e-6+1.74251e-6im
 -0.000171078-1.62315e-20im     4.90318e-6+1.51933e-6im
 -0.000172917+6.92626e-20im     4.80709e-6+1.34041e-6im
 -0.000174289-1.93188e-20im  …   4.7313e-6+1.20609e-6im
 -0.000175201+1.07167e-20im     4.67897e-6+1.11652e-6im
 -0.000175656+2.03995e-20im     4.65226e-6+1.07173e-6im

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
```jldoctest
julia> using ContinuousWavelets, Random

julia> rng = MersenneTwister(23425); Y = randn(rng, 2053, 4);

julia> X = Y .+ 3;

julia> c = wavelet(morl)

julia> wCo = waveletCoherence(X, Y, c);
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.038173051785201154
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6

julia> wCo[:,:,1,1]
2053×18 Matrix{Float64}:
 0.969189   0.854792  0.974984  …  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.96828    0.853079  0.974938     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.966475   0.849751  0.974852     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.963798   0.844997  0.974734     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.960287   0.839082  0.974599     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.955994   0.832323  0.974462  …  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.950986   0.825068  0.974338     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.945346   0.817671  0.974243     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.939174   0.810471  0.974189     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.932583   0.80377   0.974186     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 ⋮                              ⋱                      ⋮
 0.162165   0.991927  0.904406     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.12791    0.992445  0.900202  …  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0984489  0.99293   0.896122     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.074116   0.99337   0.892289     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0549532  0.993756  0.888831     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0407007  0.994076  0.885876     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0308526  0.994324  0.883542  …  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0247819  0.994494  0.881927     1.0  1.0  1.0  1.0  1.0  1.0  1.0
 0.0219171  0.994579  0.881101     1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> wCo[:,:,1,2]
2053×18 Matrix{Float64}:
 0.952459  0.57212   0.844804  0.64797   …  0.199122  0.0702893   0.325306
 0.951494  0.57209   0.844596  0.645967     0.19763   0.0657064   0.320263
 0.949561  0.572297  0.844265  0.641983     0.195414  0.0571731   0.310414
 0.946653  0.573206  0.843967  0.636057     0.193827  0.0458668   0.29622
 0.942758  0.575366  0.843905  0.628223     0.194496  0.0333503   0.278331
 0.937858  0.579285  0.844304  0.618499  …  0.198969  0.0213685   0.257544
 0.931926  0.585335  0.845377  0.606863     0.20843   0.0116364   0.234747
 0.924926  0.593703  0.847297  0.593255     0.22351   0.00564711  0.210866
 0.916807  0.604378  0.850178  0.577569     0.244245  0.00452665  0.186805
 0.907506  0.617179  0.854061  0.559668     0.270138  0.00894809  0.16339
 ⋮                                       ⋱  ⋮
 0.826733  0.989244  0.326962  0.841909     0.563908  0.285219    0.289197
 0.849596  0.989932  0.351184  0.853911  …  0.561157  0.284495    0.290488
 0.869187  0.990576  0.373794  0.864666     0.559331  0.285274    0.291523
 0.885661  0.991161  0.394437  0.874143     0.558378  0.287236    0.292392
 0.899158  0.991673  0.412703  0.882275     0.558148  0.289962    0.293152
 0.909802  0.992099  0.428133  0.888969     0.558419  0.292967    0.293819
 0.917694  0.992429  0.44025   0.894117  …  0.558932  0.295761    0.294374
 0.922912  0.992654  0.448617  0.897618     0.55944   0.297905    0.294779
 0.925507  0.992768  0.452893  0.89939      0.559747  0.299066    0.294994

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
