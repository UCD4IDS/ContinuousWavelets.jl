function morsefreq(c::CWT{W,T,Morse,N}) where {W,T,N}
    
    # measures of frequency for generalized Morse wavelet. [with F. Rekibi]
    # the output returns the modal or peak
    
    # For be=0, the "wavelet" becomes an analytic lowpass filter
    
    
    # Lilly and Olhede (2009).  Higher-order properties of analytic wavelets.  
    # IEEE Trans. Sig. Proc., 57 (1), 146--160.


    ga = c.waveType.ga;
    be = c.waveType.be;
    
    fm = exp.((log.(be) - log.(ga)) ./ ga);
    
    if sum(be.==0) != 0 && size(fm) == ()
        fm = (log(2))^(1 / ga); 
    elseif sum(be.==0) != 0 && size(fm) != ()
        fm[be.==0] = (log(2))^(1 / ga[be.==0]); 
    end

    fm = fm / (2 * pi);
    
    return fm

end


"""
    daughter = mother(this::CWT{W, T, <:ContWaveClass, N},
                        s::Real, nInOctave::Int, ω::AbstractArray{<:Real,1}) where {W, T, N}

given a CWT object, return a rescaled version of the mother wavelet, in the
fourier domain. ω is the frequency, which is fftshift-ed. s is the scale
variable.
"""
function mother(this::CWT{W,T,Morlet,N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W,T,N}
    constant = this.σ[3] * (π)^(1 / 4)
    gauss = exp.(-(this.σ[1] .- ω / s).^2 / (sWidth^2))
    shift = this.σ[2] * exp.(-1 / 2 * (ω / s).^2)
    daughter = constant .* (gauss .- shift)
    return normalize(daughter, s, this.p)
end

function mother(this::CWT{W,T,<:Paul,N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W,T,N}
    daughter = zeros(length(ω))
    constant = (2^this.α) / sqrt((this.α) * gamma(2 * (this.α)))
    polynomial = (ω[ω .>= 0] / s).^(this.α)
    expDecay = exp.(-(ω[ω .>= 0] / s))
    daughter[ω .>= 0] = constant .* polynomial .* expDecay
    return normalize(daughter, s, this.p)
end

function mother(this::CWT{W,T,<:Dog,N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W,T,N}
    constant = im^(this.α) / sqrt(gamma((this.α) + 1 / 2))
    polynomial = (ω / s).^(this.α)
    gauss = exp.(-(ω / s).^2 / 2)
    daughter =  constant .* polynomial .* gauss
    return normalize(daughter, s, this.p)
end

function mother(this::CWT{W,T,Morse,N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W,T,N}
    
    ga = this.waveType.ga;
    be = this.waveType.be;
    cf = this.waveType.cf;
    p = this.p;

    fo = morsefreq(this) * 2 * pi;
    fact = cf/fo;
    

    om = (ω / s);

    if be == 0
        daughter = 2 * exp.(-om.^ga);
    else
        daughter = 2 * exp.(-be.*log.(fo) .+ fo.^ga .+ be.*log.(om) .- om.^ga);
    end
    
    daughter[1] = 1/2*daughter[1]; # Due to unit step function
    # Ensure nice lowpass filters for beta=0;
    # Otherwise, doesn't matter since wavelets vanishes at zero frequency

    if any(daughter .!= daughter)
        @warn "the given values of gamma and beta are numerically unstable"
        daughter[daughter .!= daughter] .= 0;
    end
    
    return ContinuousWavelets.normalize(daughter, s, p)
end

function mother(this::CWT{W,T,<:ContOrtho,N}, s::Real, itpψ,
                ω::AbstractArray{<:Real,1}, n, n1) where {W,T,N}
    daughter = itpψ(range(1, stop=length(itpψ), step=length(itpψ) / n1 * s))
    daughter = padTo(daughter, n)
    daughter = circshift(daughter, -round(Int, n1 / s / 2))
    return daughter
end
function normalize(daughter, s, p)
    if p == Inf
        normTerm = maximum(abs.(daughter))
    else
        normTerm = s^(1 / p)
    end
    return daughter ./ normTerm
end

"""
    father(c::CWT, ω, averagingType::aT) where {aT<:Averaging}

this creates the averaging function, which covers the low frequency
information, and is emphatically not analytic. aT determines whether it has the
same form as the wavelets (`Father`), or just a bandpass `Dirac`.
`c.averagingLength` gives the number of octaves (base 2) that are covered by
the averaging function. The width is then derived so that it matches the next
wavelet at 1σ.

For the Morlet wavelet, the distribution is just a Gaussian, so it has variance
1/s^2 and mean σ[1]*s set the variance so that the averaging function has σ/2 at
the central frequency of the last scale.

For the Paul wavelets, it's a easy calculation to see that the mean of a paul
wavelet of order m is (m+1)/s, while σ=sqrt(m+1)/s. So we set the variance so
that the averaging function has 1σ at the central frequency of the last scale.

the derivative of a Gaussian (Dog) has a pretty nasty form for the mean and
variance; eventually, if you set 1/2 * σ_{averaging}=⟨ω⟩_{highest scale wavelet}, you
will get the scale of the averaging function to be
`s*gamma((c.α+2)/2)/sqrt(gamma((c.α+1)/2)*(gamma((c.α+3)/2)-gamma((c.α+2)/2)))`

"""
function father(c::CWT, ω, averagingType::Father, sWidth)
    s = 2^(min(c.averagingLength + getMinScaling(c) - 1, 0))
    s0, ω_shift = locationShift(c, s, ω, sWidth)
    averaging = adjust(c) .* mother(c, s0, 1, ω_shift)
end

function father(c::CWT{B,T,W}, ω, averagingType::Father,
                       fullVersion, s, N, n1) where {W <: ContOrtho,B,T}
    itp = genInterp(fullVersion)

    φ = itp(range(1, stop=length(fullVersion),
                  step=length(fullVersion) / n1 * s))
    φ = circshift(padTo(φ, N), -round(Int, n1 / s / 2))
end

# dirac version (that is, just a window around zero)
function father(c::CWT{<:WaveletBoundary,T}, ω,
                       averagingType::Dirac, sWidth) where {T}
    s = 2^(c.averagingLength) * sWidth
    averaging = zeros(T, size(ω))
    upperBound = getUpperBound(c, s)
    averaging[abs.(ω) .<= upperBound] .= 1
    return averaging
end

function father(c::CWT{W, T, <:Morse}, ω, averagingType::ContinuousWavelets.Father, sWidth) where {W,T}
    s = 2^(min(c.averagingLength + getMinScaling(c) - 1, 0)) 
    s0, ω_shift = locationShift(c, s, ω, sWidth)
    averaging = adjust(c) .* mother(c, s0, sWidth, ω_shift)
    return averaging
end



@doc """
      (daughters, ω) = computeWavelets(n1::Integer, c::CWT{W}; T=Float64, J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real,
                                                                   W<:WaveletBoundary, V}
just precomputes the wavelets used by transform c::CWT{W}. For details, see cwt
"""
function computeWavelets(n1::Integer, c::CWT{B,CT,W}; T=Float64, J1::Int64=-1, dt::S=NaN, s0::V=NaN, space=false) where {S <: Real,B <: WaveletBoundary,V,W,CT}
    # don't alter scaling with sampling information if it doesn't exists
    fλ = (4 * π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    if isnan(dt) || (dt < 0)
        dt = 1
    end
    # smallest resolvable scale
    if isnan(s0) || (s0 < 0)
        s0 = 2 * dt / fλ
    end

    nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1, c)
    # padding determines the actual number of elements
    n, nSpace = setn(n1, c)
    # indicates whether we should keep a spot for the father wavelet
    isAve = !(typeof(c.averagingType) <: NoAve)
    # I guess matlab did occasionally do something useful

    ω = computeOmega(n1, nSpace, n)
    daughters = analyticOrNot(c, n, totalWavelets)

    # if the nOctaves is small enough there are none not covered by the
    # averaging, just use that
    if totalWavelets == 1
        @warn "an averaging length of $(c.averagingLength) results in fewer octaves than the minimum frequency $(getMinScaling(c)) for this wavelet type. Either decreasing the averagingLength or increase extraOctaves"
        onlyFather = father(c, ω, c.averagingType, sWidth[1])
        return onlyFather, ω
    end

    for (curWave, s) in enumerate(sRange)
        daughters[:,curWave + isAve] = mother(c, s, sWidth[curWave], ω)
    end
    if isAve # should we include the father?
        @debug "c = $(c), $(c.averagingType), size(ω)= $(size(ω))"
        daughters[:, 1] = father(c, ω, c.averagingType, sWidth[1])# [1:(n1+1)]
    end

    # adjust by the frame bound
    if c.frameBound > 0
        daughters = daughters .* (c.frameBound / norm(daughters, 2))
    end
    testFourierDomainProperties(daughters, isAve)
    if space
        x = zeros(T, n1); x[ceil(Int, n1 / 2)] = 1
        return cwt(x, c, daughters), ω
    else
        return (daughters, ω)
    end
end

function computeWavelets(n1::Integer, c::CWT{B,CT,W}; T=Float64, space=false) where {S <: Real,B <: WaveletBoundary,V,W <: ContOrtho,CT}

    # padding determines the actual number of elements
    n, nSpace = setn(n1, c)
    nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1, c)
    # indicates whether we should keep a spot for the father wavelet
    isAve = !(typeof(c.averagingType) <: NoAve)
    # I guess matlab did occasionally do something useful


    ω = computeOmega(n1, nSpace, n)
    daughters = zeros(nSpace, totalWavelets)
    φ, ψ, ψLen = getContWaveFromOrtho(c, nSpace)
    itpψ = genInterp(ψ)
    # if the nOctaves is small enough there are none not covered by the
    # averaging, so just use that
    if totalWavelets == 1
        onlyFather = father(c, ω, c.averagingType, φ, c.averagingLength, nSpace, n1)
        return rfft(onlyFather, 1), ω
    end

    for (curWave, s) in enumerate(sRange)
        daughters[:,curWave + isAve] = mother(c, s, itpψ, ω, nSpace, n1)
    end
    if isAve
        daughters[:, 1] = father(c, ω, c.averagingType, φ, sRange[1],
                                        nSpace, n1)
    end
    # switch to fourier domain and normalize appropriately
    daughters = rfft(daughters, 1)
    repeatSRange = [sRange[1:isAve]...; sRange...]
    for i in 1:size(daughters, 2)
        daughters[:, i] = normalize(daughters[:, i],
                                            repeatSRange[i],
                                            c.p)
    end

    # adjust by the frame bound
    if c.frameBound > 0
        daughters = daughters .* (c.frameBound / norm(daughters, 2))
    end
    testFourierDomainProperties(daughters, isAve)
    if space
        x = zeros(T, n1); x[ceil(Int, n1 / 2)] = 1
        return cwt(x, c, daughters), ω
    else
        return (daughters, ω)
    end
end

# not
function analyticOrNot(c::CWT{W,T,<:Union{Dog,ContOrtho},N}, n, totalWavelets) where {W,T,N}
    # Odd derivatives (and any orthogonal wavelets) introduce an imaginary
    # term, so we need a complex representation
    if abs(c.α % 2) == 1
        daughters = zeros(Complex{T}, n, totalWavelets)
    else
        daughters = zeros(T, n, totalWavelets)
    end
    return daughters
end

function analyticOrNot(c::CWT{W,T,<:Union{Morlet,Paul,Morse},N}, n, totalWavelets) where {W,T,N}
    daughters = zeros(T, n, totalWavelets)
    return daughters
end

function computeOmega(nOriginal, nSpace, nFreq)
    range(0, nOriginal >> 1 + 1, length=nFreq) # max size is the last frequency in the rfft of the original data size
end
# convenience methods

# it's ok to just hand the total size, even if we're not transforming across
# all dimensions
function computeWavelets(Y::Tuple, c::CWT{W}; T=Float64) where {S <: Real,W <: WaveletBoundary}
    return computeWavelets(Y[1], c; T=T)
end
function computeWavelets(Y::AbstractArray{<:Integer}, c::CWT{W}; T=Float64) where {S <: Real,W <: WaveletBoundary}
    return computeWavelets(Y[1], c, T=T)
end

# also ok to just hand the whole thing being transformed
function computeWavelets(Y::AbstractArray{<:Number}, c::CWT{W}; T=Float64) where {S <: Real,W <: WaveletBoundary}
    return computeWavelets(size(Y)[1], c, T=T)
end
