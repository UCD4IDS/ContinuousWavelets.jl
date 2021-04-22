using Documenter, ContinuousWavelets
ENV["GKSwstype"] = "100"
makedocs(sitename="ContinuousWavelets.jl",
         pages=[
             "basic usage" => "index.md",
             "Install" => "installation.md",
             "CWT" => [
                 "Available Wavelet Families" => "coreType.md",
                 "CWT Construction" => "CWTConstruction.md",
                 "Wavelet Spacing" => "spacing.md",
                 "Boundary Conditions" => "bound.md",
                 "Inversion" => "inverse.md"
             ],
         ])

deploydocs(
    repo="github.com/dsweber2/ContinuousWavelets.jl.git",
)
