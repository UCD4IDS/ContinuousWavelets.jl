struct Morlet <: ContWaveClass
    σ::Float64 # \sigma is the time/space trade-off. as sigma->0, the spacial resolution increases; below 5, there is a danger of being non-analytic. Default is 5.8
    κσ::Float64
    cσ::Float64
end
"""
    morl = Morlet(σ::T) where T<: Real

    return the Morlet wavelet with first central frequency parameter σ, which controls the time-frequency trade-off. As σ goes to zero, all of the information becomes spatial. Default is `` 2π``, and it is recommended that you stay within 3-12. Above this range, the overlap between wavelets becomes impractically small. Below it, the mean subtraction term approaches the magnitude of the wavelet, so they become a sum of Gaussians rather than Gaussians.
"""
function Morlet(σ::T) where T <: Real
    κσ = exp(-σ^2 / 2)
    cσ = 1. / sqrt(1 + κσ^2 - 2 * exp(-3 * σ^2 / 4))
    Morlet(σ, κσ, cσ)
end
Morlet() = Morlet(2π)
class(::Morlet) = "Morlet"; name(::Morlet) = "morl"; vanishingmoments(::Morlet) = 0; isAnalytic(::Morlet) = true
const morl = Morlet()
Base.show(io::IO, x::Morlet) = print(io, "Morlet mean $(x.σ)")


struct Morse <: ContinuousWavelets.ContWaveClass
    ga::Float64
    be::Float64
    cf::Float64
end
"""
    morse = Morse(ga::T,be::T,cf::T) where T<: Float64

    return the Morse wavelet with the central frequency parameter cf, gamma parameter ga and beta parameter be.
"""
function Morse_convert(ga::T,be::T,cf::T) where T <: Real
    ga, be, cf = Float64.(ga), Float64.(be), Float64.(cf)
    Morse(ga,be,cf)
end

Morse() = Morse_convert(3, 1, 1)
class(::Morse) = "Morse"; name(::Morse) = "morse"; 
vanishingmoments(::Morse) = 0; isAnalytic(::Morse) = true
const morse = Morse()
Base.show(io::IO, x::Morse) = print(io, "Morse gamma = $(x.ga), Morse beta = $(x.be), center frequency = $(x.cf)")

# TODO: include a sombrero wavelet, which is dog2.

# Parameterized classes

# abstract type Dog <: ContWaveClass end
# abstract type Paul <: ContWaveClass end

# continuous parameterized
for (TYPE, NAMEBASE, MOMENTS, RANGE, ISAN) in ((:Paul, "paul", -1, 1:20, true), # moments? TODO: is this a good range of parameters?
        (:Dog, "dog",  -1, 0:6, false),)
    @eval begin
        struct $TYPE{N} <: ContWaveClass end# $TYPE end
        class(::$TYPE) = $(string(TYPE))
        name(::$TYPE{N}) where N = string($NAMEBASE, N)
        vanishingmoments(::$TYPE{N}) where N = -1
        order(::$TYPE{N}) where N = N # either order for Paul wavelets, or number of derivatives for DOGs
        Base.show(io::IO, x::$TYPE{N}) where N = print(io, $NAMEBASE * " order $(N)")
        isAnalytic(::$TYPE) = $(ISAN)
    end
    for NUM in RANGE
        CONSTNAME = Symbol(NAMEBASE, NUM)
        @eval begin
            const $CONSTNAME = $TYPE{$NUM}()                  # type shortcut
            export $CONSTNAME
    end
end
end

# this is a bit broken, but could unify these
# abstract type SupaType end
# macro makeType(FAMILY, FORPRINTING, PARAMS, RANGE)
#     parsedParams = mapreduce(x->"$(x[1])::$(x[2])", (a,b) ->a*"\n"*b, @eval($PARAMS))
#     println(parsedParams)
#     @eval begin
#         localParams =
#         struct $FAMILY{N} <: SupaType
#             :($(parsedParams))
#         end
#         family(::$FAMILY) = $(string(FAMILY))
#         order(::$FAMILY{N}) where N = N
#     end
#     # constants for convienence
#     for NUM in RANGE
#         CONSTNAME = Symbol(FAMILY, NUM)
#         @eval begin
#             const $CONSTNAME = $FAMILY{$NUM}()
#             export $CONSTNAME
#         end
#     end
# end
# adapting the orthogonal wavelet classes to the continuous case
struct ContOrtho{OWT} <: ContWaveClass where OWT <: OrthoFilter
    o::OWT
end
ContOrtho(o::W) where W <: WT.OrthoWaveletClass = ContOrtho{typeof(wavelet(o))}(wavelet(o))
class(a::ContOrtho{OWT}) where OWT = "Continuous $(class(a.o))"
name(a::ContOrtho{OWT}) where OWT = "c$(name(a.o))"
vanishingmoments(a::ContOrtho{OWT}) where OWT = vanishingmoments(a.o)
isAnalytic(::ContOrtho) = false
qmf(w::ContOrtho) = w.o.qmf
Base.show(io::IO, x::ContOrtho{W}) where W = print(io, "Continuous $(x.o.name)")

# Single classes
const cHaar = ContOrtho(WT.haar)
const cBeyl = ContOrtho(WT.beyl)
const cVaid = ContOrtho(WT.vaid)
# parametric orthogonal
for (TYPE, NAMEBASE, RANGE) in ((:Daubechies, "Db", 1:10),
        (:Coiflet, "Coif", 2:2:8),
        (:Symlet, "Sym", 4:10),
        (:Battle, "Batt", 2:2:6),)
    for NUM in RANGE
        CONSTNAME = Symbol(string("c", NAMEBASE, NUM))
        @eval begin
            const $CONSTNAME = ContOrtho(WT.$TYPE{$NUM}())        # type shortcut
            export $CONSTNAME
    end
end
end







### averaging types (for now this is already exported by wavelets, so just use that)

abstract type Average end
struct Dirac <: Average end
struct Father <: Average end
struct NoAve <: Average end

function Base.show(io::IO, cf::Dirac)
    print(io, "Dirac")
end
function Base.show(io::IO, cf::Father)
    print(io, "Father Wavelet")
end
function Base.show(io::IO, cf::NoAve)
    print(io, "No Averaging")
end
