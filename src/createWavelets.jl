"""
    daughter = Daughter(this::CFW{W, T, <:ContinuousWaveletClass, N}, 
                        s::Real, nInOctave::Int, ω::AbstractArray{<:Real,1}) where {W, T, N}

given a CFW object, return a rescaled version of the mother wavelet, in the
fourier domain. ω is the frequency, which is fftshift-ed. s is the scale 
variable.
"""
function Mother(this::CFW{W, T, Morlet, N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W, T, N}
    constant = this.σ[3]*(π)^(1/4)
    gauss = exp.(-(this.σ[1].-ω/s).^2/(sWidth^2))
    shift = this.σ[2]*exp.(-1/2*(ω/s).^2)
    daughter = constant .* (gauss .- shift) 
    return normalize(daughter, s, this.normalization)
end

function Mother(this::CFW{W, T, <:Paul, N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W, T, N}
    daughter = zeros(length(ω))
    constant = (2^this.α) / sqrt((this.α) * gamma(2*(this.α)))
    polynomial = (ω[ω.>=0]/s).^(this.α)
    expDecay = exp.(-(ω[ω.>=0]/s))
    daughter[ω.>=0]= constant .* polynomial .* expDecay
    return normalize(daughter, s, this.normalization)
end

function Mother(this::CFW{W, T, <:Dog, N}, s::Real, sWidth::Real,
                  ω::AbstractArray{<:Real,1}) where {W, T, N}
    constant = im^(this.α)/sqrt(gamma((this.α)+1/2))
    polynomial = (ω/s).^(this.α)
    gauss = exp.(-(ω/s).^2/2)
    daughter =  constant .* polynomial .* gauss
    return normalize(daughter, s, this.normalization)
end

function normalize(daughter, s, p)
    if p == Inf
        normTerm = maximum(abs.(daughter))
    else
        normTerm = s^(1/p)
    end
    return daughter ./ normTerm
end

"""
    findAveraging(c::CFW, ω, averagingType::aT) where {aT<:Averaging}

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
function findAveraging(c::CFW{W, T, WT}, ω, averagingType::Father, sWidth)
    s = 2^(c.averagingLength-1)
    s0, ω_shift = locationShift(c, s, ω, sWidth)
    averaging = adjust(c) .* Mother(c, s0, 1, ω_shift)
end

function findAveraging(c::CFW{W, T, WT}, ω, averagingType::Father, 
                       sWidth) where {WT <: ContOrtho}
    s = 2^(c.averagingLength-1)
    s0, ω_shift = locationShift(c, s, ω, sWidth)
    averaging = adjust(c) .* Mother(c, s0, 1, ω_shift)
end

# dirac version (that is, just a window around zero)
function findAveraging(c::CFW{<:WaveletBoundary, T}, ω,
                       averagingType::Dirac, sWidth) where {T}
    s = 2^(c.averagingLength) * sWidth
    averaging = zeros(T, size(ω))
    upperBound = getUpperBound(c, s)
    averaging[abs.(ω) .<= upperBound] .= 1
    return averaging
end






@doc """
      (daughters, ω) = computeWavelets(n1::Integer, c::CFW{W}; T=Float64, J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real,
                                                                   W<:WT.WaveletBoundary, V}
just precomputes the wavelets used by transform c::CFW{W}. For details, see cwt
"""
function computeWavelets(n1::Integer, c::CFW{B, CT, W}; T=Float64, J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real,
                                                                   B<:WT.WaveletBoundary, V, WT}
    # don't alter scaling with sampling information if it doesn't exists
    fλ = (4*π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    if isnan(dt) || (dt<0)
        dt = 1
    end
    # smallest resolvable scale
    if isnan(s0) || (s0<0)
        s0 = 2 * dt / fλ
    end
    # J1 is the total number of scales
    # if J1<0
    #     J1 = Int(round(log2(n1 * dt / s0) * c.scalingFactor))
    # end

    nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1,c)
    # padding determines the actual number of elements
    n = setn(n1,c)
    # indicates whether we should keep a spot for the father wavelet
    isAve = (c.averagingLength > 0 && !(typeof(c.averagingType) <: NoAve)) ? 1 : 0
    # I guess matlab did occasionally do something useful


    ω, daughters = analyticOrNot(c, n, totalWavelets)
    

    # if the nOctaves is small enough there are none not covered by the
    # averaging, just use that
    if round(nOctaves) < 0
        father = findAveraging(c,ω, c.averagingType, sWidth[1])
        return father, ω
    end

    for (curWave, s) in enumerate(sRange)
        daughters[:,curWave+isAve] = Mother(c, s, sWidth[curWave], ω)
    end
    if isAve>0 # should we include the father?
        @debug "c = $(c), $(c.averagingType), size(ω)= $(size(ω))"
        daughters[:, 1] = findAveraging(c, ω, c.averagingType, sWidth[1])#[1:(n1+1)]
    end

    # adjust by the frame bound
    if c.frameBound > 0
        daughters = daughters.*(c.frameBound/norm(daughters, 2))
    end
    testFourierDomainProperties(daughters, isAve)
    return (daughters, ω)
end

function computeWavelets(n1::Integer, c::CFW{B, CT, W}; T=Float64, J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real,
                                                                   B<:WT.WaveletBoundary, V, WT<:ContOrtho}
    # don't alter scaling with sampling information if it doesn't exists
    fλ = (4*π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    if isnan(dt) || (dt<0)
        dt = 1
    end
    # smallest resolvable scale
    if isnan(s0) || (s0<0)
        s0 = 2 * dt / fλ
    end

    # padding determines the actual number of elements
    n = setn(n1,c)
    nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1, c)
    # indicates whether we should keep a spot for the father wavelet
    isAve = (c.averagingLength > 0 && !(typeof(c.averagingType) <: NoAve)) ? 1 : 0
    # I guess matlab did occasionally do something useful


    ω, daughters = analyticOrNot(c, n, totalWavelets)
    

    # if the nOctaves is small enough there are none not covered by the
    # averaging, just use that
    if round(nOctaves) < 0
        father = findAveraging(c,ω, c.averagingType, sWidth[1])
        return father, ω
    end

    for (curWave, s) in enumerate(sRange)
        daughters[:,curWave+isAve] = Mother(c, s, sWidth[curWave], ω)
    end
    if isAve>0 # should we include the father?
        @debug "c = $(c), $(c.averagingType), size(ω)= $(size(ω))"
        daughters[:, 1] = findAveraging(c, ω, c.averagingType, sWidth[1])#[1:(n1+1)]
    end

    # adjust by the frame bound
    if c.frameBound > 0
        daughters = daughters.*(c.frameBound/norm(daughters, 2))
    end
    testFourierDomainProperties(daughters, isAve)
    return (daughters, ω)
end
# not
function analyticOrNot(c::CFW{W, T, <:Dog, N}, n, totalWavelets) where {W,T,N}
    ω = (0:(n-1))*2π
    # Odd derivatives introduce an imaginary term, so we need a complex representation
    if c.α % 2 == 1
        daughters = zeros(Complex{T}, n, totalWavelets)
    else
        daughters = zeros(T, n, totalWavelets)
    end
    return (ω, daughters)
end

function analyticOrNot(c::CFW{W, T, <:Union{Morlet, Paul}, N}, n, totalWavelets) where {W,T,N}
    ω = (0:(n-1))*2π
    daughters = zeros(T, n, totalWavelets)
    return (ω, daughters)
end



# convenience methods

# it's ok to just hand the total size, even if we're not transforming across
# all dimensions
function computeWavelets(Y::Tuple, c::CFW{W}; T=Float64) where {S<:Real,
                                                     W<:WT.WaveletBoundary}
    return computeWavelets(Y[1], c; T=T)
end
function computeWavelets(Y::AbstractArray{<:Integer}, c::CFW{W}; T=Float64) where {S<:Real,
                                                     W<:WT.WaveletBoundary}
    return computeWavelets(Y[1], c, T=T)
end

# also ok to just hand the whole thing being transformed
function computeWavelets(Y::AbstractArray{<:Number}, c::CFW{W}; T=Float64) where {S<:Real,
                                                     W<:WT.WaveletBoundary}
    return computeWavelets(size(Y)[1], c, T=T)
end

