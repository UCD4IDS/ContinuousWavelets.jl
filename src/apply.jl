# CWT (continuous wavelet transform)
# cwt(Y::AbstractVector, ::ContinuousWavelet)

@doc """
     wave = cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WT}, daughters, rfftPlan =
             plan_rfft([1]), fftPlan = plan_fft([1])) where {N, T<:Real,
                                                             S<:Real,
                                                             U<:Number,
                                                             W<:WT.WaveletBoundary,
                                                             WT<:Union{<:WT.Morlet,
                                                                       <:WT.Paul}}

  return the continuous wavelet transform along the first axis with averaging.
  `wave`, is (signalLength)×(nscales+1)×(previous dimensions), of type T of
  Y. averagingLength defines the number of octaves (powers of 2) that are
  replaced by an averaging function. This has form averagingType, which can be
  one of `Mother()` or `Dirac()`- in the `Mother()` case, it uses the same form
  as for the wavelets, while the `Dirac` uses a constant window. J1 is the
  total number of scales; default (when J1=NaN, or is negative) is just under
  the maximum possible number, i.e. the log base 2 of the length of the signal,
  times the number of wavelets per octave. If you have sampling information,
  you will need to scale wave by δt^(1/2).

  """
function cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WaTy}, daughters, rfftPlan::AbstractFFTs.Plan =
             plan_rfft([1]), fftPlan = plan_fft([1])) where {N, T<:Real,
                                                             S<:Real,
                                                             W<:WT.WaveletBoundary,
                                                             WaTy<:Union{<:WT.Morlet,
                                                                         <:WT.Paul}}
    # This is for analytic wavelets, so we need to treat the positive and
    # negative frequencies differently, even for real data

    # TODO: complex input version of this
    @assert typeof(N)<:Integer
    # vectors behave a bit strangely, so we reshape them
    if N==1
        Y= reshape(Y,(length(Y), 1))
    end
    n1 = size(Y, 1);
    
    _, nScales, _ = getNWavelets(n1, c)
    @debug "" nScales
    #....construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(c)()) #this function is defined below

    # check if the plans we were given are dummies or not
    if size(rfftPlan)==(1,)
        rfftPlan = plan_rfft(x, 1)
    end
    if size(fftPlan)==(1,)
        fftPlan = plan_fft(x, 1)
    end
    n = size(x, 1)
    
    x̂ = rfftPlan * x
    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging
    if nScales <= 0 || size(daughters,2) == 1
        daughters = daughters[:,1:1]
        nScales = 1
    end

    wave = zeros(Complex{T}, size(x)..., nScales);  # result array;
    # faster if we put the example index on the outside loop through all scales
    # and compute transform
    actuallyTransform!(wave, daughters,x̂, fftPlan, c.waveType, c.averagingType)
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave)-1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...] 
    if N==1
        wave = dropdims(wave, dims=3)
    end

    return wave
end


function cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WaTy}, daughters, rfftPlan =
             plan_rfft([1])) where {N, T<:Real, S<:Real, U<:Number,
                                    W<:WT.WaveletBoundary, WaTy<:WT.Dog}
    # Dog doesn't need a fft because it is strictly real
    # TODO: complex input version of this
    @assert typeof(N)<:Integer
    # vectors behave a bit strangely, so we reshape them
    if N==1
        Y= reshape(Y,(length(Y), 1))
    end

    n1 = size(Y, 1);
    
    _, nScales, _ = getNWavelets(n1, c)
    #....construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(c)())

    # check if the plans we were given are dummies or not
    if size(rfftPlan)==(1,)
        rfftPlan = plan_rfft(x, 1)
    end
    n = size(x, 1)

    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging. Or if there's only averaging
    if nScales <= 1 || size(daughters,2) == 1
        daughters = daughters[:,1:1]
        nScales = 1
    end

    x̂ = rfftPlan * x
    
    wave = zeros(Complex{T}, size(x)..., nScales);  # result array;
    # faster if we put the example index on the outside loop through all scales
    # and compute transform


    actuallyTransform!(wave, daughters,x̂, rfftPlan, c.waveType)
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave)-1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...] 

    if N==1
        wave = dropdims(wave, dims=3)
    end

    return real.(wave)
end




function reflect(Y, bt)
    n1 = size(Y, 1)
    if bt == WT.padded
        base2 = round(Int,log(n1)/log(2));   # power of 2 nearest to N
        x = cat(Y, zeros(2^(base2+1)-n1, size(Y)[2:end]...), dims=1)
    elseif bt == WT.DEFAULT_BOUNDARY
        x = cat(Y, reverse(Y,dims = 1), dims = 1)
    else
        x = Y
    end
    return x
end

function actuallyTransform!(wave, daughters, x̂, fftPlan, analytic::Union{<:WT.Morlet,
                                                                         <:WT.Paul},
                            averagingType::Union{WT.Father, WT.Dirac})
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    isSourceOdd = mod(size(wave,1)+1,2)
    # the averaging function isn't analytic, so we need to do both positive and
    # negative frequencies
    tmpWave = x̂ .* daughters[:,1]
    wave[(n1+1):end, outer..., 1] = reverse(conj.(tmpWave[2:end-isSourceOdd,
                                                          outer...]),dims=1)
    wave[1:n1, outer..., 1] = tmpWave
    wave[:, outer..., 1] = fftPlan \ (wave[:, outer..., 1])  # wavelet transform
    for j in 2:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end
function actuallyTransform!(wave, daughters, x̂, fftPlan, analytic::Union{<:WT.Morlet,
                                                                         <:WT.Paul},
                            averagingType::WT.NoAve)
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    # the no averaging version
    for j in 1:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

function actuallyTransform!(wave, daughters, x̂, rfftPlan, analytic::Union{<:WT.Dog})
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    for j in 1:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = rfftPlan \ (wave[1:n1, outer..., j])  # wavelet transform
    end
end



function cwt(Y::AbstractArray{T}, c::CFW{W}; varArgs...) where {T<:Number, S<:Real, V<: Real,
                                                                                        W<:WT.WaveletBoundary}
    daughters,ω = computeWavelets(size(Y, 1), c; varArgs...) 
    return cwt(Y, c, daughters)
end


"""
period,scale, coi = caveats(Y::AbstractArray{T}, c::CFW{W}; J1::S=NaN) where {T<:Real, S<:Real, W<:WT.WaveletBoundary}

returns the period, the scales, and the cone of influence for the given wavelet transform. If you have sampling information, you will need to scale the vector scale appropriately by 1/δt, and the actual transform by δt^(1/2).
"""
function caveats(n1, c::CFW{W}; J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real, W<:WT.WaveletBoundary, V <: Real}
    # don't alter scaling with sampling information if it doesn't exists
    fλ = (4*π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    if isnan(dt) || (dt<0)
        dt = 1
    end

    # smallest resolvable scale
    if isnan(s0) || (s0<0)
        s0 = 2 * dt / fλ
    end
    sj = s0 * 2.0.^(collect(0:J1)./c.scalingFactor)
    # Fourier equivalent frequencies
    freqs = 1 ./ (fλ .* sj)

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with
    # non-zero end-points.
    coi = (n1 / 2 .- abs.(collect(0:n1-1) .- (n1 - 1) ./ 2))
    coi = (fλ * dt / sqrt(2)).*coi


    n1 = length(Y);
    # J1 is the total number of elements
    if isnan(J1) || (J1<0)
        J1=floor(Int,(log2(n1))*c.scalingFactor);
    end
    #....construct time series to analyze, pad if necessary
    if boundaryType(c) == WT.ZPBoundary
        base2 = round(Int,log(n1)/log(2));   # power of 2 nearest to N
        n = length(Y)+2^(base2+1)-n1
    elseif boundaryType(c) == WT.PerBoundary
        n = length(Y)*2
    end
    ω = [0:floor(Int, n/2); -floor(Int,n/2)+1:-1]*2π
    period = c.fourierFactor*2 .^((0:J1)/c.scalingFactor)
    scale = [1E-5; 1:((n1+1)/2-1); reverse((1:(n1/2-1)),dims=1); 1E-5]
    coi = c.coi*scale  # COI [Sec.3g]
    return sj, freqs, period, scale, coi
end
cwt(Y::AbstractArray{T}, w::WT.ContinuousWaveletClass; J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {T<:Real, S<:Real, V<:Real} = cwt(Y,CFW(w),J1=J1,dt=dt,s0=s0)
caveats(Y::AbstractArray{T}, w::WT.ContinuousWaveletClass; J1::S=NaN) where {T<: Real, S<: Real} = caveats(Y,CFW(w),J1=J1)
cwt(Y::AbstractArray{T}) where T<:Real = cwt(Y,WT.Morlet())
caveats(Y::AbstractArray{T}) where T<:Real = caveats(Y,WT.Morlet())

"""
    icwt(W::AbstractArray{T}, c::CFW{W}, sj::AbstractArray; dt::S=NaN,
         dj::V=1/12) where {S<:Real, V<:Real, W<:WT.WaveletBoundary}
    icwt(W::AbstractArray{T}, c::CFW{W}; dt::S=NaN, s0::S
         dj::V=1/12) where {T<:Real, S<:Real, V<:Real,
                            W<:WT.WaveletBoundary} 


return the inverse continuous wavelet transform
"""
function icwt(W::AbstractArray, c::CFW, sj::AbstractArray; dt::S=NaN,
              dj::V=1/12) where {S<:Real, V<:Real} 
    if isnan(dt) || (dt<0)
        dt = 1
    end

    # Torrence and Compo (1998), eq. (11)
    n = WT.setn(size(W,1), c)
    ω = (0:(n-1))*2π
    ψ = WT.Mother(c, 1, 1, ω)[1:(end-1)]
    println("size(sj) = $(size(sj))")
    println("size(W) = $(size(W))")
    println("size(ψ) = $(size(ψ))")
    iW = (dj * sqrt(dt) / 0.776 * psi(c,0)) .* sum((real.(W) ./ sqrt.(sj')), dims=2)

    return iW
end
function icwt(W::AbstractArray, c::CFW; dt::Real=NaN, dj::Real=1/12, J1::Real=NaN)

    fλ = (4*π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    n1 = size(W, 1);
    # J1 is the total number of elements
    if isnan(J1) || (J1<0)
        J1=floor(Int,(log2(n1))*c.scalingFactor);
    end

    sj =  WT.getScales(n1, c)

    if isnan(dt) || (dt<0)
        dt = 1
    end

    return icwt(W, c, sj, dt=dt, dj=1/(c.scalingFactor))
end
icwt(Y::AbstractArray, w::WT.ContinuousWaveletClass; dj::T=1/12, dt::S=NaN,
     s0::V=NaN) where {S<:Real, T<:Real, V<:Real} = icwt(Y,CFW(w))
icwt(Y::AbstractArray) = icwt(Y,WT.Morlet())

function psi(c::CFW{W}, t::Int64) where W<:WT.WaveletBoundary
    return real.(π^(-0.25) * exp.(im*c.σ[1]*t - t^2 / 2))
end




# CWT (continuous wavelet transform directly) TODO: direct if sufficiently small

# TODO: continuous inverse, when defined
#icwt(::AbstractVector, ::ContinuousWavelet)
