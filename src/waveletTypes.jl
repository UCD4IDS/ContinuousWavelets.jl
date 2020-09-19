


struct Morlet <: ContinuousWaveletClass
    σ::Float64 # \sigma is the time/space trade-off. as sigma->0, the spacial resolution increases; below 5, there is a danger of being non-analytic. Default is 5.8
    κσ::Float64
    cσ::Float64
end
"""
    morl = Morlet(σ::T) where T<: Real

    return the Morlet wavelet with first central frequency parameter σ, which controls the time-frequency trade-off. As σ goes to zero, all of the information becomes spatial. Default is `` 2π``, and it is recommended that you stay within 3-12. Above this range, the overlap between wavelets becomes impractically small. Below it, the mean subtraction term approaches the magnitude of the wavelet, so they become a sum of Gaussians rather than Gaussians.
"""
function Morlet(σ::T) where T<:Real
    κσ=exp(-σ^2/2)
    cσ=1. /sqrt(1+κσ^2-2*exp(-3*σ^2/4))
    Morlet(σ,κσ,cσ)
end
Morlet() = Morlet(2π)
class(::Morlet) = "Morlet"; name(::Morlet) = "morl"; vanishingmoments(::Morlet)=0
const morl = Morlet()
Base.show(io::IO, x::Morlet) = print(io,"Morlet mean $(x.σ)")

# TODO: include a sombrero wavelet, which is dog2.

# Parameterized classes

# abstract type Dog <: ContinuousWaveletClass end
# abstract type Paul <: ContinuousWaveletClass end

# continuous parameterized
for (TYPE, NAMEBASE, MOMENTS, RANGE) in (
        (:Paul, "paul", -1, 1:20), # moments? TODO: is this a good range of parameters?
        (:Dog, "dog",  -1, 0:6), # moments?
        )
    @eval begin
        struct $TYPE{N} <: ContinuousWaveletClass end#$TYPE end
        class(::$TYPE) = $(string(TYPE))
        name(::$TYPE{N}) where N = string($NAMEBASE,N)
        vanishingmoments(::$TYPE{N}) where N = -1
        order(::$TYPE{N}) where N = N # either order for Paul wavelets, or number of derivatives for DOGs
        Base.show(io::IO, x::$TYPE{N}) where N = print(io,$NAMEBASE * " order $(N)")
    end
    for NUM in RANGE
        CONSTNAME = Symbol(NAMEBASE,NUM)
        @eval begin
            const $CONSTNAME = $TYPE{$NUM}()                  # type shortcut
        end
    end
end




# adapting the orthogonal wavelet classes to the continuous case
struct ContOrtho{OWT} <: ContinuousWaveletClass
    o::OWT
end
ContOrtho(o::WT) where WT <: Union{OrthoWaveletClass,BiOrthoWaveletClass} = ContOrtho{WT}(o)
class(a::ContOrtho{OWT}) = "Continuous $(class(a.o))"
name(a::ContOrtho{OWT}) = "c$(name(a.o))"
vanishingmoments(a::ContOrtho{OWT}) = vanishingmoments(a.o)
# Single classes
const cHaar = ContOrtho(WT.haar)
const cBeyl = ContOrtho(WT.beyl)
const cVaid = ContOrtho(WT.vaid)

for (TYPE, NAMEBASE, RANGE) in (
        (:Daubechies, "Db", 1:10),
        (:Coiflet, "Coif", 2:2:8),
        (:Symlet, "Sym", 4:10),
        (:Battle, "Batt", 2:2:6),
        )
    for NUM in RANGE
        CONSTNAME = Symbol(string("c", NAMEBASE, NUM))
        @eval begin
            const $CONSTNAME = ContOrtho(WT.$TYPE{$NUM}())        # type shortcut
        end
    end
end







### averaging types

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
