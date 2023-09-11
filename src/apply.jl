# CWT (continuous wavelet transform)
# cwt(Y::AbstractVector, ::ContWave)

@doc """
     wave = cwt(Y::AbstractArray{T,N}, c::CWT{W, S, WT}, daughters, rfftPlan =
             plan_rfft([1]), fftPlan = plan_fft([1])) where {N, T<:Real,
                                                             S<:Real,
                                                             U<:Number,
                                                             W<:WaveletBoundary,
                                                             WT<:ContWaveClass}

  return the continuous wavelet transform along the first axis with averaging.
  `wave`, is (signalLength)×(nscales+1)×(previous dimensions), of type `T` of
  `Y`. `averagingLength` defines the number of octaves (powers of 2) that are
  replaced by an averaging function. This has form `averagingType`, which can be
  one of `Father()` or `Dirac()`- in the `Father()` case, it uses the same form
  as for the wavelets, while the `Dirac` uses a constant window. If you have
  sampling information, you will need to scale wave by δt^(1/p). The default
  assumption is that the sampling rate is 2kHz.

"""
function cwt(Y::AbstractArray{T,N}, cWav::CWT, daughters, fftPlans = 1) where {N,T}
    @assert typeof(N) <: Integer
    # vectors behave a bit strangely, so we reshape them
    if N == 1
        Y = reshape(Y, (length(Y), 1))
    end
    n1 = size(Y, 1)

    _, nScales, _ = getNWavelets(n1, cWav)
    #construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(cWav)()) #this function is defined below

    # check if the plans we were given are dummies or not
    x̂, fromPlan = prepSignalAndPlans(x, cWav, fftPlans)
    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging
    if nScales <= 0 || size(daughters, 2) == 1
        daughters = daughters[:, 1:1]
        nScales = 1
    end

    if isAnalytic(cWav.waveType)
        OutType = ensureComplex(T)
    else
        OutType = T
    end

    wave = zeros(OutType, size(x)..., nScales)  # result array
    # faster if we put the example index on the outside loop through all scales
    # and compute transform
    if isAnalytic(cWav.waveType)
        if eltype(x) <: Real
            analyticTransformReal!(wave, daughters, x̂, fromPlan, cWav.averagingType)
        else
            analyticTransformComplex!(wave, daughters, x̂, fromPlan, cWav.averagingType)
        end
    else
        otherwiseTransform!(wave, daughters, x̂, fromPlan, cWav.averagingType)
    end
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave)-1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...]
    if N == 1
        wave = dropdims(wave, dims = 3)
    end

    return wave
end

function ensureComplex(T)
    if T <: Real
        return Complex{T}
    else
        return T
    end
end

# there are 4 cases to deal with
#       input Type | Real | Complex
#       analytic?  |----------------
#             yes  | both | fft
#              no  | rfft | fft
#  Analytic on Real input
function prepSignalAndPlans(x::AbstractArray{T},
    cWav::CWT{W,S,WaTy,N,true},
    fftPlans) where {T<:Real,W,S,WaTy,N}
    # analytic wavelets that are being applied on real inputs
    if fftPlans isa Tuple{<:AbstractFFTs.Plan{<:Real},<:AbstractFFTs.Plan{<:Complex}}
        # they handed us the right kind of thing, so no need to make new ones
        x̂ = fftPlans[1] * x
        fromPlan = fftPlans[2]
    else
        toPlan = plan_rfft(x, 1)
        x̂ = toPlan * x
        fromPlan = plan_fft(x, 1)
    end
    return x̂, fromPlan
end

#  Non-analytic on Real input
function prepSignalAndPlans(x::AbstractArray{T},
    cWav::CWT{W,S,WaTy,N,false},
    fftPlans) where {T<:Real,W,S,WaTy,N}
    # real wavelets that are being applied on real inputs
    if fftPlans isa AbstractFFTs.Plan{<:Real}
        # they handed us the right kind of thing, so no need to make new ones
        x̂ = fftPlans * x
        fromPlan = fftPlans
    else
        fromPlan = plan_rfft(x, 1)
        x̂ = fromPlan * x
    end
    return x̂, fromPlan
end
# complex input
function prepSignalAndPlans(x::AbstractArray{T}, cWav, fftPlans) where {T<:Complex}
    # real wavelets that are being applied on real inputs
    if fftPlans isa AbstractFFTs.Plan{<:Complex}
        # they handed us the right kind of thing, so no need to make new ones
        x̂ = fftPlans * x
        fromPlan = fftPlans
    else
        fromPlan = plan_fft(x, 1)
        x̂ = fromPlan * x
    end
    return x̂, fromPlan
end

# analytic on real data with an averaging function
function analyticTransformReal!(wave, daughters, x̂, fftPlan, ::Union{Father,Dirac})
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    isSourceEven = mod(size(wave, 1) + 1, 2)
    # the averaging function isn't analytic, so we need to do both positive and
    # negative frequencies
    @views tmpWave = x̂ .* daughters[:, 1]
    @views wave[(n1+1):end, outer..., 1] = reverse(conj.(tmpWave[2:end-isSourceEven,
            outer...]),
        dims = 1)
    @views wave[1:n1, outer..., 1] = tmpWave
    @views wave[:, outer..., 1] = fftPlan \ (wave[:, outer..., 1])  # averaging
    for j = 2:size(daughters, 2)
        @views wave[1:n1, outer..., j] = x̂ .* daughters[:, j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

# analytic on complex data with an averaging function
function analyticTransformComplex!(wave, daughters, x̂, fftPlan, ::Union{Father,Dirac})
    outer = axes(x̂)[2:end]
    n1 = size(daughters, 1)
    isSourceEven = mod(size(wave, 1) + 1, 2)
    # the averaging function isn't analytic, so we need to do both positive and
    # negative frequencies
    @views positiveFreqs = x̂[1:n1, outer...] .* daughters[:, 1]
    @views negativeFreqs = x̂[(n1-isSourceEven+1):end, outer...] .*
                           reverse(conj.(daughters[2:end, 1]))
    @views wave[(n1-isSourceEven+1):end, outer..., 1] = negativeFreqs
    @views wave[1:n1, outer..., 1] = positiveFreqs
    @views wave[:, outer..., 1] = fftPlan \ (wave[:, outer..., 1])  # averaging
    for j = 2:size(daughters, 2)
        @views wave[1:n1, outer..., j] = x̂[1:n1, outer...] .* daughters[:, j]
        @views wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

function analyticTransformComplex!(wave, daughters, x̂, fftPlan, averagingType)
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    for j = 1:size(daughters, 2)
        @views wave[1:n1, outer..., j] = x̂[1:n1, outer...] .* daughters[:, j]
        @views wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

# analytic on real data without an averaging function
function analyticTransformReal!(wave, daughters, x̂, fftPlan, ::NoAve)
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    # the no averaging version
    for j = 1:size(daughters, 2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:, j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end


function otherwiseTransform!(wave::AbstractArray{<:Real},
    daughters,
    x̂,
    fromPlan,
    averagingType)
    # real wavelets on real data: that just makes sense
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    for j = 1:size(daughters, 2)
        @views tmp = x̂ .* daughters[:, j]
        @views wave[:, outer..., j] = fromPlan \ tmp  # wavelet transform
    end
end

# if it isn't analytic, the output is complex only if the input is complex
function otherwiseTransform!(wave::AbstractArray{<:Complex},
    daughters,
    x̂,
    fromPlan,
    averagingType)
    # applying a real transform to complex data is maybe a bit odd, but you do you
    outer = axes(x̂)[2:end]
    n1 = size(daughters, 1)
    isSourceEven = mod(size(fromPlan, 1) + 1, 2)
    for j = 1:size(daughters, 2)
        @views wave[1:n1, outer..., j] = @views x̂[1:n1, outer...] .* daughters[:, j]
        @views wave[n1-isSourceEven+1:end, outer..., j] = x̂[n1-isSourceEven+1:end,
            outer...] .* reverse(conj.(daughters[2:end,
            j]))
        @views wave[:, outer..., j] = fromPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

function reflect(Y, bt)
    n1 = size(Y, 1)
    if typeof(bt) <: ZPBoundary
        base2 = ceil(Int, log2(n1))   # power of 2 nearest to N
        x = cat(Y, zeros(2^(base2) - n1, size(Y)[2:end]...), dims = 1)
    elseif typeof(bt) <: SymBoundary
        x = cat(Y, reverse(Y, dims = 1), dims = 1)
    else
        x = Y
    end
    return x
end


function cwt(Y::AbstractArray{T},
    c::CWT{W};
    varArgs...) where {T<:Number,W<:WaveletBoundary}
    daughters, ω = computeWavelets(size(Y, 1), c; varArgs...)
    return cwt(Y, c, daughters)
end

function cwt(Y::AbstractArray{T}, w::ContWaveClass; varargs...) where {T<:Number}
    cwt(Y, CWT(w); varargs...)
end
cwt(Y::AbstractArray{T}) where {T<:Real} = cwt(Y, Morlet())

abstract type InverseType end
struct DualFrames <: InverseType end
struct NaiveDelta <: InverseType end
struct PenroseDelta <: InverseType end

"""
    icwt(res::AbstractArray, cWav::CWT, inverseStyle::InverseType=PenroseDelta())
Compute the inverse wavelet transform using one of three dual frames. The default uses delta functions with weights chosen via a least squares method, the `PenroseDelta()` below. This is chosen as a default because the Morlet wavelets tend to fail catastrophically using the canonical dual frame (the `dualFrames()` type).

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::PenroseDelta)
Return the inverse continuous wavelet transform, computed using the simple dual frame ``β_jδ_{ji}``, where ``β_j`` is chosen to solve the least squares problem ``\\|Ŵβ-1\\|_2^2``, where ``Ŵ`` is the Fourier domain representation of the `cWav` wavelets. In both this case and `NaiveDelta()`, the fourier transform of ``δ`` is the constant function, thus this least squares problem.

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::NaiveDelta)
Return the inverse continuous wavelet transform, computed using the simple dual frame ``β_jδ_{ji}``, where ``β_j`` is chosen to negate the scale factor ``(^1/_s)^{^1/_p}``. Generally less accurate than choosing the weights using `PenroseDelta`. This is the method discussed in Torrence and Compo.

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::dualFrames)
Return the inverse continuous wavelet transform, computed using the canonical dual frame ``\\tilde{\\widehat{ψ}} = \\frac{ψ̂_n(ω)}{∑_n\\|ψ̂_n(ω)\\|^2}``. The algorithm is to compute the cwt again, but using the canonical dual frame; consequentially, it is the most computationally intensive of the three algorithms, and typically the best behaved. Will be numerically unstable if the high frequencies of all of the wavelets are too small however, and tends to fail spectacularly in this case.

"""
function icwt(res::AbstractArray, cWav::CWT, ::PenroseDelta)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    β = computeDualWeights(Ŵ, cWav)
    testDualCoverage(β, Ŵ)
    compXRecon = sum(res .* β, dims = 2)
    imagXRecon = irfft(im * rfft(imag.(compXRecon), 1), size(compXRecon, 1)) # turns out the dual frame for the imaginary part is rather gross in the time domain
    return imagXRecon + real.(compXRecon)
end

function icwt(res::AbstractArray, cWav::CWT, ::NaiveDelta)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    β = computeNaiveDualWeights(Ŵ, cWav, size(res, 1))
    testDualCoverage(β, Ŵ)
    compXRecon = sum(res .* β, dims = 2)
    imagXRecon = irfft(im * rfft(imag.(compXRecon), 1), size(compXRecon, 1)) # turns out the dual frame for the imaginary part is rather gross in the time domain
    return imagXRecon + real.(compXRecon)
end

function icwt(res::AbstractArray, cWav::CWT, ::DualFrames)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    canonDualFrames = explicitConstruction(Ŵ)
    testDualCoverage(canonDualFrames)
    convolved = cwt(res, cWav, canonDualFrames)
    ax = axes(convolved)
    @views xRecon = sum(convolved[:, i, i, ax[4:end]...] for i = 1:size(Ŵ, 2))
    return xRecon
end

function icwt(Y::AbstractArray, w::ContWaveClass; varargs...)
    icwt(Y, CWT(w), PenroseDelta(); varargs...)
end
icwt(Y::AbstractArray; varargs...) = icwt(Y, Morlet(), PenroseDelta(); varargs...)

# CWT (continuous wavelet transform directly) TODO: direct if sufficiently small
