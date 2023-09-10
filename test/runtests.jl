using ContinuousWavelets, Wavelets, Interpolations, LinearAlgebra
using Test, Documenter
using FFTW
using Logging, Random
inGithubAction = get(() -> "", ENV, "JULIA_IN_GITHUB_ACTION") == "true"
inGithubActionOnMac = get(() -> "", ENV, "JULIA_IN_GITHUB_ACTION_ON_MAC") == "macOS-latest"
# these make sure that the printing width/length is kept to a reasonable amout for actually reading the docs
ENV["LINES"] = "9"
ENV["COLUMNS"] = "60"
@testset "ContinuousWavelets.jl" begin
    if (inGithubAction && !inGithubActionOnMac)
        doctest(ContinuousWavelets)
    end
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
using LanguageServer
using JuliaFormatter
format(".")
