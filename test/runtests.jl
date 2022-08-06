using ContinuousWavelets, Wavelets, Interpolations, LinearAlgebra
using Test, Documenter
using FFTW
using Logging, Random
@testset "ContinuousWavelets.jl" begin
    doctest(ContinuousWavelets)
    include("basicTypesAndNumber.jl")
    include("deltaSpikes.jl")
    include("utilsTests.jl")
    include("defaultProperties.jl")
    include("inversionTests.jl")
end
# TODO:
#       test averaging types
#            various extra dimensions
#            inverse is actually functional
