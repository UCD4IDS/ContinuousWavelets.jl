using ContinuousWavelets, Wavelets, Interpolations, LinearAlgebra
using Test, Documenter
using FFTW
using Logging, Random
ENV["LINES"] = "9"
ENV["COLUMNS"] = "60"
@testset "ContinuousWavelets.jl" begin
    doctest(ContinuousWavelets, doctestfilters=[
        r"\@ ContinuousWavelets .*",
        r"[ +-][0-9]\.[0-9]{3,5}e-1[5-9]",
        r"[ +-][0-9]\.[0-9]{3,5}e-[2-9][0-9]",
        r"im {2,7}",
    ])
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
