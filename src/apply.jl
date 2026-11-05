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
    _checkMatchingDevice(Y, daughters)
    # vectors behave a bit strangely, so we reshape them
    if N == 1
        Y = reshape(Y, (length(Y), 1))
    end
    n1 = size(Y, 1)

    _, nScales, _ = getNWavelets(n1, cWav)
    #construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(cWav)()) #this function is defined below

    # check if the plans we were given are dummies or not. Only the forward
    # plan is used here; the inverse is built against the batched shape below.
    x̂, _ = prepSignalAndPlans(x, cWav, fftPlans)
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

    # One plan for the entire wavelet bank rather than one FFT per wavelet
    invPlan = prepBatchedInversePlan(x, cWav, nScales, fftPlans)

    # Allocate directly in the final layout, at the final (untrimmed) length 
    wave = similar(x, OutType, (n1, nScales, size(x)[2:end]...))

    if isAnalytic(cWav.waveType)
        if eltype(x) <: Real
            analyticTransformReal!(wave, daughters, x̂, invPlan, cWav.averagingType)
        else
            analyticTransformComplex!(wave, daughters, x̂, invPlan, cWav.averagingType)
        end
    else
        otherwiseTransform!(wave, daughters, x̂, invPlan, cWav.averagingType)
    end

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


_batchedSize(x::AbstractArray, nScales) = (size(x, 1), nScales, size(x)[2:end]...)

_isPlanOfSize(p, sz) = (p isa AbstractFFTs.Plan) && (size(p) == sz)

# analytic wavelets on real input: complex plan, used in the inverse direction
function prepBatchedInversePlan(x::AbstractArray{T},
    cWav::CWT{W,S,WaTy,N,true},
    nScales, fftPlans) where {T<:Real,W,S,WaTy,N}
    sz = _batchedSize(x, nScales)
    # a caller passing (rfftPlan, fftPlan) may supply an already-batched second
    # element and skip construction entirely; a signal-shaped one is ignored
    if fftPlans isa Tuple{<:AbstractFFTs.Plan{<:Real},<:AbstractFFTs.Plan{<:Complex}} &&
       _isPlanOfSize(fftPlans[2], sz)
        return fftPlans[2]
    end
    return plan_fft(similar(x, complex(float(T)), sz), 1)
end

# non-analytic wavelets on real input
function prepBatchedInversePlan(x::AbstractArray{T},
    cWav::CWT{W,S,WaTy,N,false},
    nScales, fftPlans) where {T<:Real,W,S,WaTy,N}
    return plan_rfft(similar(x, float(T), _batchedSize(x, nScales)), 1)
end

# complex input, either flavour of wavelet
function prepBatchedInversePlan(x::AbstractArray{T}, cWav, nScales,
    fftPlans) where {T<:Complex}
    return plan_fft(similar(x, complex(float(real(T))), _batchedSize(x, nScales)), 1)
end


function _freqBuffer(x̂::AbstractArray, sz, fullyWritten::Bool)
    buf = similar(x̂, sz)
    fullyWritten || fill!(buf, 0)
    return buf
end

# insert a singleton scale axis as dimension 2, so that a (nFreq × nScales)
# daughter bank broadcasts across it
_withScaleAxis(x̂::AbstractArray) = reshape(x̂, size(x̂, 1), 1, size(x̂)[2:end]...)

# copy the leading (unpadded) rows of the batched transform into the result
function _trimTo!(wave::AbstractArray, padded::AbstractArray)
    @views wave .= padded[1:size(wave, 1), ntuple(_ -> Colon(), ndims(wave) - 1)...]
    return wave
end


# analytic on real data with an averaging function
function analyticTransformReal!(wave, daughters, x̂, fftPlan, ::Union{Father,Dirac})
    outer = axes(x̂)[2:end]
    nFreq = size(x̂, 1)
    nPad = size(fftPlan, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    isSourceEven = mod(nPad + 1, 2)
    negFreqEnd = nFreq - isSourceEven
    X = _withScaleAxis(x̂)

    Ẑ = _freqBuffer(x̂, (nPad, nScales, size(x̂)[2:end]...), false)
    # positive frequencies, every wavelet at once
    @views Ẑ[1:nFreq, 1:nD, outer...] .= X .* daughters
    # the averaging function isn't analytic, so it -- and only it -- also needs
    # the negative frequencies. Source rows 2:negFreqEnd and destination rows
    # nFreq+1:end are disjoint, so reading and writing Ẑ here is safe.
    @views Ẑ[(nFreq+1):end, 1, outer...] .= conj.(Ẑ[negFreqEnd:-1:2, 1, outer...])

    _trimTo!(wave, fftPlan \ Ẑ)
    return wave
end

# analytic on real data without an averaging function
function analyticTransformReal!(wave, daughters, x̂, fftPlan, ::NoAve)
    outer = axes(x̂)[2:end]
    nFreq = size(x̂, 1)
    nPad = size(fftPlan, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    X = _withScaleAxis(x̂)

    Ẑ = _freqBuffer(x̂, (nPad, nScales, size(x̂)[2:end]...), false)
    @views Ẑ[1:nFreq, 1:nD, outer...] .= X .* daughters

    _trimTo!(wave, fftPlan \ Ẑ)
    return wave
end

# analytic on complex data with an averaging function
function analyticTransformComplex!(wave, daughters, x̂, fftPlan, ::Union{Father,Dirac})
    outer = axes(x̂)[2:end]
    nFreq = size(daughters, 1)
    nPad = size(fftPlan, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    isSourceEven = mod(nPad + 1, 2)
    negFreqStart = nFreq - isSourceEven + 1
    X = _withScaleAxis(x̂)

    Ẑ = _freqBuffer(x̂, (nPad, nScales, size(x̂)[2:end]...), false)
    @views Ẑ[1:nFreq, 1:nD, outer...] .= X[1:nFreq, :, outer...] .* daughters
    @views Ẑ[negFreqStart:end, 1, outer...] .= X[negFreqStart:end, 1, outer...] .*
                                               conj.(daughters[nFreq:-1:2, 1])

    _trimTo!(wave, fftPlan \ Ẑ)
    return wave
end

# analytic on complex data without an averaging function
function analyticTransformComplex!(wave, daughters, x̂, fftPlan, averagingType)
    outer = axes(x̂)[2:end]
    nFreq = size(daughters, 1)
    nPad = size(fftPlan, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    X = _withScaleAxis(x̂)

    Ẑ = _freqBuffer(x̂, (nPad, nScales, size(x̂)[2:end]...), false)
    @views Ẑ[1:nFreq, 1:nD, outer...] .= X[1:nFreq, :, outer...] .* daughters

    _trimTo!(wave, fftPlan \ Ẑ)
    return wave
end

function otherwiseTransform!(wave::AbstractArray{<:Real},
    daughters,
    x̂,
    fromPlan,
    averagingType)
    # real wavelets on real data: that just makes sense
    outer = axes(x̂)[2:end]
    nFreq = size(x̂, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    X = _withScaleAxis(x̂)

    # rfft domain: nFreq rows in, the plan's inverse produces the nPad rows
    Ẑ = _freqBuffer(x̂, (nFreq, nScales, size(x̂)[2:end]...), nD == nScales)
    @views Ẑ[:, 1:nD, outer...] .= X .* daughters

    _trimTo!(wave, fromPlan \ Ẑ)
    return wave
end

# if it isn't analytic, the output is complex only if the input is complex
function otherwiseTransform!(wave::AbstractArray{<:Complex},
    daughters,
    x̂,
    fromPlan,
    averagingType)
    # applying a real transform to complex data is maybe a bit odd, but you do you
    outer = axes(x̂)[2:end]
    nFreq = size(daughters, 1)
    nPad = size(fromPlan, 1)
    nScales = size(wave, 2)
    nD = size(daughters, 2)
    isSourceEven = mod(nPad + 1, 2)
    negStart = nFreq - isSourceEven + 1
    X = _withScaleAxis(x̂)

    # rows 1:nFreq and negStart:nPad together cover every row, so the buffer
    # only needs zeroing when there are more scales than daughters
    Ẑ = _freqBuffer(x̂, (nPad, nScales, size(x̂)[2:end]...), nD == nScales)
    @views Ẑ[1:nFreq, 1:nD, outer...] .= X[1:nFreq, :, outer...] .* daughters
    @views Ẑ[negStart:end, 1:nD, outer...] .= X[negStart:end, :, outer...] .*
                                              conj.(daughters[nFreq:-1:2, 1:nD])

    _trimTo!(wave, fromPlan \ Ẑ)
    return wave
end

function reflect(Y, bt)
    n1 = size(Y, 1)
    if typeof(bt) <: ZPBoundary
        base2 = ceil(Int, log2(n1))   # power of 2 nearest to N
        x = cat(Y, zeros(eltype(Y), 2^(base2) - n1, size(Y)[2:end]...), dims = 1)
    elseif typeof(bt) <: SymBoundary
        x = cat(Y, _reverseDim1(Y), dims = 1)
    else
        x = Y
    end
    return x
end

_matchPrecision(reference::AbstractArray, x::AbstractArray) =
    eltype(x) <: Complex ? Complex{real(eltype(reference))}.(x) : real(eltype(reference)).(x)

_matchDevice(reference::AbstractArray, x::AbstractArray) = adapt(_backend(reference), x)

function cwt(Y::AbstractArray{T}, c::CWT{W}; varArgs...) where {T<:Number,W<:WaveletBoundary}
    daughters, ω = computeWavelets(size(Y, 1), c; varArgs...)
    daughters = _matchPrecision(Y, daughters)
    daughters = _matchDevice(Y, daughters)
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

#= Remedy for an error that emerged from GPU implementation. _backend was not 
   returning the compute device as expected. To remedy, recurses on `parent`
   covers for every wrapper that defines it — Base's views, reshapes and reinterprets,
   LinearAlgebra's Adjoint/Transpose, and package types such as OffsetArrays —
   and leaves anything that does not define `parent` reporting its own type. =#
function _unwrapped(x::AbstractArray)
    p = parent(x)
    return p === x ? x : _unwrapped(p)
end

_backend(x::AbstractArray) = Base.typename(typeof(_unwrapped(x))).wrapper

function _checkMatchingDevice(Y::AbstractArray, daughters::AbstractArray)
    _backend(Y) === _backend(daughters) || error(
        "cwt: signal is on $(_backend(Y)) but daughters is on " *
        "$(_backend(daughters)). Put both on the same device: move the signal " *
        "with Array(Y) or CuArray(Y), or pass daughters yourself via " *
        "cwt(Y, cWav, daughters).")
    return nothing
end

_reverseDim1(A::AbstractVector) = A[length(A):-1:1]
_reverseDim1(A::AbstractArray) = A[size(A, 1):-1:1, ntuple(_ -> Colon(), ndims(A) - 1)...]

"""
    icwt(res::AbstractArray, cWav::CWT, inverseStyle::InverseType=PenroseDelta())
Compute the inverse wavelet transform using one of three dual frames. The default uses delta functions with weights chosen via a least squares method, the `PenroseDelta()` below. This is chosen as a default because the Morlet wavelets tend to fail catastrophically using the canonical dual frame (the `dualFrames()` type).

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::PenroseDelta)
Return the inverse continuous wavelet transform, computed using the simple dual frame ``β_jδ_{ji}``, where ``β_j`` is chosen to solve the least squares problem ``\\|Ŵβ-1\\|_2^2``, where ``Ŵ`` is the Fourier domain representation of the `cWav` wavelets. In both this case and `NaiveDelta()`, the fourier transform of ``δ`` is the constant function, thus this least squares problem.

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::NaiveDelta)
Return the inverse continuous wavelet transform, computed using the simple dual frame ``β_jδ_{ji}``, where ``β_j`` is chosen to negate the scale factor ``(^1/_s)^{^1/_p}``. Generally less accurate than choosing the weights using `PenroseDelta`. This is the method discussed in Torrence and Compo.

    icwt(res::AbstractArray, cWav::CWT, inverseStyle::DualFrames)
Return the inverse continuous wavelet transform, computed using the canonical dual frame ``\\tilde{\\widehat{ψ}} = \\frac{ψ̂_n(ω)}{∑_n\\|ψ̂_n(ω)\\|^2}``. The algorithm is to compute the cwt again, but using the canonical dual frame; consequentially, it is the most computationally intensive of the three algorithms, and typically the best behaved. Will be numerically unstable if the high frequencies of all of the wavelets are too small however, and tends to fail spectacularly in this case.

"""
function icwt(res::AbstractArray, cWav::CWT, ::PenroseDelta)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    Ŵ = _matchPrecision(res, Ŵ)
    β = computeDualWeights(Ŵ, cWav)
    β = _matchPrecision(res, β)
    testDualCoverage(β, Ŵ)
    β = _matchDevice(res, β)
    compXRecon = sum(res .* β, dims = 2)
    imagXRecon = irfft(im * rfft(imag.(compXRecon), 1), size(compXRecon, 1))
    return imagXRecon + real.(compXRecon)
end

function icwt(res::AbstractArray, cWav::CWT, ::NaiveDelta)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    Ŵ = _matchPrecision(res, Ŵ)
    β = computeNaiveDualWeights(Ŵ, cWav, size(res, 1))
    β = _matchPrecision(res, β)
    testDualCoverage(β, Ŵ)
    β = _matchDevice(res, β)
    compXRecon = sum(res .* β, dims = 2)
    imagXRecon = irfft(im * rfft(imag.(compXRecon), 1), size(compXRecon, 1))
    return imagXRecon + real.(compXRecon)
end

function icwt(res::AbstractArray, cWav::CWT, ::DualFrames)
    Ŵ = computeWavelets(size(res, 1), cWav)[1]
    Ŵ = _matchPrecision(res, Ŵ)
    canonDualFrames = explicitConstruction(Ŵ)
    canonDualFrames = _matchPrecision(res, canonDualFrames)
    testDualCoverage(canonDualFrames)
    canonDualFrames = _matchDevice(res, canonDualFrames)
    convolved = cwt(res, cWav, canonDualFrames)
    ax = axes(convolved)
    @views xRecon = sum(convolved[:, i, i, ax[4:end]...] for i = 1:size(Ŵ, 2))
    return xRecon
end

function icwt(Y::AbstractArray, w::ContWaveClass; varargs...)
    icwt(Y, CWT(w), PenroseDelta(); varargs...)
end
icwt(Y::AbstractArray; varargs...) = icwt(Y, Morlet(), PenroseDelta(); varargs...)

# CWT (continuous wavelet transform directly) TODO: direct if sufficiently small