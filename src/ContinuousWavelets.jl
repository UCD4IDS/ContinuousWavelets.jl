module ContinuousWavelets
using Wavelets
using Interpolations
using AbstractFFTs
using FFTW
using LinearAlgebra

import Wavelets.wavelet, Wavelets.getNWavelets, Wavelets.computeWavelets # yo ho yo ho
import Wavelets.WT.name, Wavelets.WT.name, Wavelets.WT.class,
   Wavelets.WT.vanishingmoments 

export ContWave 

# Boundaries
export WaveletBoundary, PerBoundary, ZPBoundary, NullBoundary, SymBoundary,
    Periodic, DEFAULT_BOUNDARY, padded, NaivePer, SymBound
# waveletTypes (note there are also export statements in the for loops)
export Morlet, Paul, Dog, ContOrtho, name, vanishingmoments, ContOrtho
export cHaar, cBeyl, cVaid
# averaging types
export Average, Dirac, Father, NoAve
# CWT constructors
export wavelet, waveletType
# general utilities
export qmf, computeWavelets


abstract type ContWave{Boundary,T} end # equivalent to ContinuousWavelet in Wavelets.jl
# BOUNDARY TYPES

abstract type WaveletBoundary end
# periodic (default)
struct PerBoundary <: WaveletBoundary end
# zero padding
struct ZPBoundary <: WaveletBoundary end
# constant padding
#struct CPBoundary <: WaveletBoundary end
struct NullBoundary <: WaveletBoundary end
# symmetric boundary (as in the DCTII)
struct SymBoundary <: WaveletBoundary end
# and so on...


const Periodic = PerBoundary()
const DEFAULT_BOUNDARY = PerBoundary()
const padded = ZPBoundary()
const NaivePer = NullBoundary()
const SymBound = SymBoundary()

abstract type ContWaveClass end # equivalent to ContinuousWaveletClass in Wavelets.jl
include("waveletTypes.jl")
include("CWTConstruction.jl")
include("utils.jl")
include("sanityChecks.jl")
include("createWavelets.jl")
include("apply.jl")
end
