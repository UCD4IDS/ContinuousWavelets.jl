using Revise
using ContinuousWavelets, Wavelets, Interpolations, LinearAlgebra
using Test
using FFTW
using Logging
@testset "ContinuousWavelets.jl" begin
    # Write your tests here.
    include("basicTypesAndNumber.jl")
    include("deltaSpikes.jl")
    include("defaultProperties.jl")
    include("inversionTests.jl")
end
# TODO:
#       test averaging types
#            various extra dimensions
#            inverse is actually functional
sj, freqs, period, scale, coi = ContinuousWavelets.caveats(128, wavelet(morl))
plot(coi)
plot(randn(10))
typeof(wavelet(morl))
