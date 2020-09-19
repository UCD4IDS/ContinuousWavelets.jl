module ContinuousWavelets
using Wavelets
using Interpolations

abstract type ContinuousWavelet{Boundary,T} end
# Write your package code here.

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

abstract type ContinuousWaveletClass end
include("waveletTypes.jl")
include("CFWConstruction.jl")
include("utils.jl")
include("createWavelets.jl")
include("apply.jl")
include("sanityChecks.jl")
end
