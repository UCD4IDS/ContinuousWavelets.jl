using Documenter, ContinuousWavelets

makedocs(sitename="ContinuousWavelets.jl",
         pages = [
             "basic usage" => "index.md",
             "Install" => "installation.md",
             "CWT" => [
                 "Available Wavelets" => "coreType.md",
                 "CWT Type" => "CWTConstruction.md",
                 "Wavelet Spacing" => "spacing.md",
                 "Boundary Conditions" => "bound.md",
             ],
         ])

deploydocs(
    repo = "github.com/dsweber2/FourierFilterFlux.jl.git",
)
