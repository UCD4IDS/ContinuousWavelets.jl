struct CWT{B, S, W<:ContWaveClass, N} <: ContWave{B, S}
    scalingFactor::S # the number of wavelets per octave, ie the scaling
                           # is s=2^(j/scalingfactor)
    decreasing::S # the amount that scalingFactor decreases per octave,
                        # to a minimum of 1
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
    normalization  ::N       # the normalization that is preserved for the
                             # wavelets as the scale changes. The conjugate
                             # p-norm is preserved for the  signal. Should be
                             # larger than 1
end

""" 
    CWT(wave::WC, scalingFactor::S=8.0, averagingType::Symbol=:Father,
        boundary::T=DEFAULT_BOUNDARY, averagingLength::Int =
        floor(Int,2*scalingFactor), frameBound::Float64=1.0,
        normalization::Float=Inf) where {WC<:ContWaveClass,
        T<:WaveletBoundary, S<:Real}

The constructor for the general type. `w` is a type of continuous wavelet,
scalingFactor is the number of wavelets between the octaves ``2^J`` and
``2^{J+1}`` (defaults to 8, which is most appropriate for music and other
audio). As this leads to excessively many high scale wavelets, `decreasing`
gives the amount that scalingFactor decreases per octave. The default boundary
condition is `periodic`, which is implemented by appending a flipped version of
the vector at the end (to eliminate edge discontinuities). Alternatives are
`ZPBoundary`, which pads with enough zeros to get to the nearest power of 2
(here the results returned by caveats are relevant, see Torrence and Compo
'97), and `NullBoundary`, which assumes the data is inherently periodic.

`averagingLength` and `averagingType` determine how wide scale information is
accounted for. `averagingLength` gives the number of wavelet octaves that are
covered by the averaging, while averaging type determines whether it is a
window `Dirac` or a wavelet `Father`. `frameBound` gives the total norm of the
whole collection, corresponding to the upper frame bound. `normalization`
refers to which p-norm is preserved as the scale changes. `normalization==2` is
the default scaling, while `normalization==Inf` gives all the same maximum
value, thus acting more like windows.
"""
function CWT(wave::WC, scalingFactor=8, boundary::T=DEFAULT_BOUNDARY,
             averagingType::A = NoAve(),
             averagingLength::Int = 0,
             frameBound=-1, normalization::N=Inf, 
             decreasing=1) where {WC<:ContWaveClass, A <: Average,
                                  T <: WaveletBoundary, N <: Real}
    @assert scalingFactor > 0
    @assert normalization >= 1
    nameWavelet = name(wave)[1:3]
    tdef = calculateProperties(wave)

    if averagingLength <= 0 || typeof(averagingType) <: NoAve
        averagingLength = 0
        averagingType = NoAve()
    end

    # S is the most permissive type of the listed variables
    S = promote_type(typeof(scalingFactor), typeof(decreasing),
                     typeof(frameBound), typeof(normalization),
                     typeof(tdef[1]), typeof(tdef[2]))
    return CWT{T, S, WC, N}(S(scalingFactor), S(decreasing), tdef...,
                            averagingLength, averagingType, S(frameBound),
                            S(normalization))
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

function eltypes(::CWT{W, T, WT, N}) where {W, T, WT, N}
    T
end
function boundaryType(::CWT{W, T, WT, N}) where {W, T, WT, N}
    W
end
function waveletType(::CWT{W, T, WT, N}) where {W, T, WT, N}
    WT
end

function Base.show(io::IO, cf::CWT{W,S,WT,N}) where {W,S,WT,N}
    print("CWT[$(cf.waveType), $(cf.averagingType), decreasing rate = "*
          "$(cf.decreasing), aveLen = $(cf.averagingLength), frame = "*
          "$(cf.frameBound), norm=$(cf.normalization)]")
end




function wavelet(cw::T; s::S=8.0, boundary::WaveletBoundary=DEFAULT_BOUNDARY,
                 averagingType::A=Father(), averagingLength::Int = 4,
                 frameBound=1, normalization::N=Inf,
                 decreasing=4) where {T<:ContWaveClass, A<:Average,
                                      S<:Real, N<:Real} 
    if typeof(s) <: AbstractFloat
        decreasing = S(decreasing)
    end
    if typeof(decreasing) <: AbstractFloat
        s = typeof(decreasing)(s)
    end

    return CWT(cw,s, boundary, averagingType, averagingLength, frameBound,
               normalization, decreasing)
end

function wavelet(c::T, boundary::WaveletBoundary) where T<:ContWaveClass
    CWT(c,8, boundary)
end
