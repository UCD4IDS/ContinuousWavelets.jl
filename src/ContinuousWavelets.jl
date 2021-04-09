module ContinuousWavelets
using Wavelets
using Interpolations
using AbstractFFTs
using FFTW
using LinearAlgebra
using SpecialFunctions

import Wavelets.wavelet, Wavelets.computeWavelets, Wavelets.cwt # yo ho yo ho
import Wavelets.WT.name, Wavelets.WT.name, Wavelets.WT.class,
   Wavelets.WT.vanishingmoments

export ContWave, CWT, cwt, icwt

# Boundaries
export WaveletBoundary, PerBoundary, ZPBoundary, NullBoundary, SymBoundary,
    Periodic, DEFAULT_BOUNDARY, padded, NaivePer, SymBound
# waveletTypes (note there are also export statements in the for loops of waveletTypes.jl)
export Morlet, Paul, Dog, ContOrtho, name, vanishingmoments, ContOrtho
export cHaar, cBeyl, cVaid, morl
# averaging types
export Average, Dirac, Father, NoAve
# CWT constructors
export wavelet, waveletType
# general utilities
export qmf, computeWavelets, getNWavelets, mother, father, getMeanFreq

"""
    ContWave{Boundary,T}
The abstract type encompassing the various types of wavelets implemented in
the package. The abstract type has parameters `Boundary<:WaveletBoundary` and
`T<:Number`, which gives the element output type. Each has both a constructor,
and a default predefined entry. These
are:
- `Morlet`: A complex approximately analytic wavelet that is just a frequency
    domain Gaussian with mean subtracted. `Morlet(σ::T) where T<: Real`. `σ`
    gives the frequency domain variance of the mother Wavelet. As `σ` goes to
    zero, all of the information becomes spatial. Default is `morl` which has
    σ=2π.
- `Paul{N}`: A complex analytic wavelet. `pauln` for n in `1:20` e.g. `paul5`
- `Dog{N}`: Derivative of a Gaussian, where N is the number of
    derivatives. `dogn` for `n` in `0:6`. The Sombrero/mexican hat/Marr wavelet
    is `n=2`.
- `ContOrtho{OWT}`. OWT is some orthogonal wavelet of type `OrthoWaveletClass`
    from [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl). This uses an
    explicit construction of the mother wavelet for these orthogonal wavelets
    to do a continuous transform. Constructed via `ContOrtho(o::W)` where `o`
    is from Wavelets.jl. Alternatively, you can get them directly as
    `ContOrtho` objects via:
   + `cHaar` Haar Wavelets
   + `cBeyl` Beylkin Wavelets
   + `cVaid` Vaidyanthan Wavelets
   + `cDbn` Daubhechies Wavelets. n ranges from `1:Inf`
   + `cCoifn` Coiflets. n ranges from `2:2:8`
   + `cSymn` Symlets. n ranges from `4:10`
   + `cBattn` Battle-Lemarie wavelets. n ranges from `2:2:6`
"""
abstract type ContWave{Boundary,T} end # equivalent to ContinuousWavelet in Wavelets.jl
# BOUNDARY TYPES

"""
the abstract type for the various types of boundaries
"""
abstract type WaveletBoundary end
# periodic (default)
"""
    PerBoundary() <: WaveletBoundary
standard periodic boundary assumption made by the fft. Aliases of `NaivePer`
and `Periodic`.
"""
struct PerBoundary <: WaveletBoundary end
# zero padding
"""
    ZPBoundary() <: WaveletBoundary
zero pads the signal before doing an fft, rounding up to the nearest power of
two. Alias of `padded`.
"""
struct ZPBoundary <: WaveletBoundary end
# constant padding
# struct CPBoundary <: WaveletBoundary end
# symmetric boundary (as in the DCTII)
"""
    SymBoundary() <: WaveletBoundary
symmetric boundary, as in the DCT type II. Repeats the edge value, which
alleviates derivative discontinuities. Alias of `DEFAULT_BOUNDARY` and
`SymBound`.
"""
struct SymBoundary <: WaveletBoundary end
# and so on...


const Periodic = PerBoundary()
const DEFAULT_BOUNDARY = SymBoundary()
const padded = ZPBoundary()
const NaivePer = PerBoundary()
const SymBound = SymBoundary()

abstract type ContWaveClass end # equivalent to ContinuousWaveletClass in Wavelets.jl
include("waveletTypes.jl")
include("CWTConstruction.jl")
include("utils.jl")
include("sanityChecks.jl")
include("createWavelets.jl")
include("apply.jl")
end
