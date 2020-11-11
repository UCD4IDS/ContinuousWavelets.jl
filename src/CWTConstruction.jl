struct CWT{B, S, W<:ContWaveClass, N} <: ContWave{B, S}
    Q::S # the number of wavelets per octave, ie the scaling
                           # is s=2^(j/scalingfactor)
    β::S # the amount that Q decreases per octave
    fourierFactor::S
    coi          ::S
    α            ::Int   # the order for a Paul and the number of derivatives
                         # for a DOG
    σ            ::Array{S} # the morlet wavelet parameters
                            # (σ,κσ,cσ). NaN if not morlet.
    waveType     ::W        # because multiple dispatch is good TODO this is actually redundant, make a function called waveType that extracts this from above
    averagingLength::Int # the number of scales to override with averaging. If
                         # you want no averaging set it to zero
    averagingType  ::Average # either Dirac or mother; the first uniformly
                             # represents the lowest frequency information,
                             # while the second constructs 
                             # a wavelet using Daughter that has mean frequency
                             # zero, and std equal to 
                             # the first non-removed wavelet's mean
    frameBound     ::S       # if positive, set the frame bound of the
                             # transform to be frameBound. Otherwise leave it
                             # so that each wavelet has an L2 norm of 1 
    p::N                     # the normalization that is preserved for the
                             # wavelets as the scale changes. The conjugate
                             # p-norm is preserved for the  signal. Should be
                             # larger than 1
end

# aliased = ((:Q,:s,:scalingFactor), (:β,:decreasing), (:p, :normalization))
function processKeywordArgs(Q,β,p;kwargs...)
    keysVarg = keys(kwargs)
    if :s in keysVarg
        Q = kwargs[:s]
        if :scalingFactor in keysVarg
            error("you used two names for Q")
        end
    elseif :scalingFactor in keysVarg
        Q = kwargs[:scalingFactor]
    end

    if :decreasing in keysVarg
        β = kwargs[:decreasing]
    end

    if :normalization in keysVarg
        p = kwargs[:normalization]
    end
    return Q,β,p
end


@doc """
    CWT(wave::WC, Q=8, boundary::T=SymBoundary(),
        averagingType::A = Father(),
        averagingLength::Int = 4, frameBound=1, p::N=Inf,
        β=4) where {WC<:ContWaveClass, A <: Average,
                    T <: WaveletBoundary, N <: Real}
"""
function CWT(wave::WC, Q=8, boundary::T=DEFAULT_BOUNDARY,
             averagingType::A = Father(),
             averagingLength::Int = 4,
             frameBound=1, p::N=Inf, 
             β=4, kwargs...) where {WC<:ContWaveClass, A <: Average,
                                  T <: WaveletBoundary, N <: Real}
    Q,β,p = processKeywordArgs(Q, β, p; kwargs...) # some names are redundant
    @assert β > 0
    @assert p >= 1
    nameWavelet = name(wave)[1:3]
    tdef = calculateProperties(wave)

    if averagingLength <= 0 || typeof(averagingType) <: NoAve
        averagingLength = 0
        averagingType = NoAve()
    end

    # S is the most permissive type of the listed variables
    S = promote_type(typeof(Q), typeof(β),
                     typeof(frameBound), typeof(p),
                     typeof(tdef[1]), typeof(tdef[2]))
    return CWT{T, S, WC, N}(S(Q), S(β), tdef...,
                            averagingLength, averagingType, S(frameBound),
                            S(p))
end
function calculateProperties(w::Morlet)
    σ = w.σ
    fourierFactor = (4*π)/(σ + sqrt(2 + σ.^2))
    coi = 1 / sqrt(2)
    α = -1
    σ = [w.σ, w.κσ, w.cσ]
    return fourierFactor, coi, α, σ, w
end

function calculateProperties(w::Dog)
    α = order(w)
    fourierFactor = 2*π*sqrt(2 ./(2*α+1))
    coi = 1/sqrt(2)
    σ = [NaN]
    return fourierFactor, coi, α, σ, w
end

function calculateProperties(w::Paul)
    α = order(w)
    fourierFactor = 4*π/(2*α+1)
    coi = sqrt(2)
    σ = [NaN]
    return fourierFactor, coi, α, σ, w
end

function calculateProperties(w::ContOrtho)
    α = -1
    fourierFactor = NaN
    coi = NaN
    σ = [NaN]
    return fourierFactor, coi, α, σ, w
end

name(s::CWT) = name(s.waveType)

function eltypes(::CWT{B, T,W, N}) where {B, T,W, N}
    T
end
function boundaryType(::CWT{B, T,W, N}) where {B, T,W, N}
    B
end
function waveletType(::CWT{B, T,W, N}) where {B, T,W, N}
    W
end

function Base.show(io::IO, cf::CWT{W,S,WT,N}) where {W,S,WT,N}
    print("CWT{$(cf.waveType), $(cf.averagingType), Q=$(Q), β=$(cf.β),"*
          "aveLen=$(cf.averagingLength), frame="* "$(cf.frameBound), norm=$(cf.p)}")
end


"""
wavelet(wave::WC, Q=8, boundary::T=DEFAULT_BOUNDARY,
                 averagingType::A = Father(), averagingLength::Int = 4,
                 frameBound=1, p::N=Inf, β=4,
                 kwargs...) where {WC<:ContWaveClass, A <: Average,
                                   T <: WaveletBoundary, N <: Real}
A constructor for the CWT type, using keyword rather than positional options.
"""
function wavelet(wave::WC; Q=8, boundary::T=DEFAULT_BOUNDARY,
                 averagingType::A = Father(), averagingLength::Int = 4,
                 frameBound=1, p::N=Inf, β=4,
                 kwargs...) where {WC<:ContWaveClass, A <: Average,
                                   T <: WaveletBoundary, N <: Real}
    return CWT(wave, Q, boundary, averagingType, averagingLength, frameBound,
               p, β, kwargs...)
end
