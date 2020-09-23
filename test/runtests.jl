using ContinuousWavelets, Wavelets, Interpolations, LinearAlgebra
using Test
using FFTW
using Logging
@testset "ContinuousWavelets.jl" begin
    # Write your tests here.
    include("basicTypesAndNumber.jl")
    include("deltaSpikes.jl")
end
# TODO: test actual values, e.g. delta spike versus generating the wavelets
#       test averaging types
#            various extra dimensions
#            inverse is actually functional
