```@meta
DocTestFilters = [r"\@ ContinuousWavelets .*", r"[ +-][0-9]\.[0-9]{5}e-[0-9][5-9]"]
```

# ContinuousWavelets

[![Build Status](https://travis-ci.com/dsweber2/ContinuousWavelets.jl.svg?branch=master)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
[![Coverage](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dsweber2.github.io/ContinuousWavelets.jl/dev/)

This package is an offshoot of [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) for the continuous wavelets.
Thanks to [Felix Gerick](https://github.com/fgerick) for the initial implementation there, with extension and further adaptation by David Weber and any other contributors listed on the right.
Currently, it implements 1D continuous wavelet transforms with the following mother wavelets:

![Mothers](docs/mothers.svg)

Which covers several standard continuous wavelet families, both real and analytic, as well as continuous versions of the orthogonal wavelet transforms implemented in [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl).

## Basic Usage

Install via the package manager and load with `using`

```julia
julia> Pkg.add("ContinuousWavelets")
julia> using ContinuousWavelets
```

Basic usage example on a doppler test function.

```jldoctest ex
julia> using Random

julia> Random.seed!(1234);

julia> using ContinuousWavelets, Wavelets

julia> n = 2047;

julia> t = range(0, n / 1000, length=n); # 1kHz sampling rate

julia> f = testfunction(n, "Doppler");

julia> c = wavelet(Morlet(π), β=2)

julia> res = ContinuousWavelets.cwt(f, c)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.061863
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:6
2047×31 Matrix{ComplexF64}:
 -2.84943e-6+3.86436e-19im  …  0.000109884+9.67013e-5im
 -2.84699e-6-6.11361e-20im     -8.24222e-5+0.000130545im
            ⋮               ⋱             ⋮
 0.000165633+5.67188e-19im  …  -2.46986e-6-1.96881e-8im
 0.000165423+1.25225e-18im     -2.63276e-6+4.61939e-8im


```

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c)

# output

2047×27 Matrix{Float64}:
 0.00262067   -0.00150482  …  -3.11967e-8  -2.84659e-8
 0.00262075   -0.00150528     -2.27191e-8  -1.96643e-8
 ⋮                         ⋱   ⋮
 0.000848453  -2.53208e-6  …  -3.11482e-9  -2.69659e-9
 0.000848413  -2.52757e-6     -4.30053e-9  -3.92791e-9
```

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
dropped = ContinuousWavelets.icwt(res, c, DualFrames())
round.(dropped,sigdigits=12)

# output

┌ Warning: the canonical dual frame is off by 3.8053445963886525e6, consider using one of the delta dual frames
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:37
2047-element Vector{Float64}:
  0.0213534682635
  0.0213551067073
  ⋮
 -0.00614525648078
 -0.00614518997054
```

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```jldoctest ex
julia> using Wavelets

julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0)

julia> res = circshift(ContinuousWavelets.cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.26778
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:10
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
 1.89367e-5  0.000266033  …  4.6727e-5    2.99983e-6
 8.33321e-5  0.000266913     1.56557e-5  -4.46419e-5
 ⋮                        ⋱  ⋮
 2.24677e-6  0.00198709   …  4.24042e-6   3.80685e-6
 2.63848e-6  0.00198004      4.3791e-6    3.47575e-6

[:, :, 2] =
  7.81007e-18  0.0226754  0.00955729  …   3.68809e-18
 -3.47114e-18  0.022684   0.00950443     -3.47114e-18
  ⋮                                   ⋱
 -9.29595e-18  0.0341512  0.0108039   …  -3.84208e-19
  1.27592e-18  0.0342157  0.0107729      -1.7043e-18

[:, :, 3] =
 -4.2736e-7   0.0059687   …  4.47839e-8  1.86209e-8
 -4.39691e-7  0.00596762     3.30771e-8  7.78201e-9
  ⋮                       ⋱  ⋮
 -9.41123e-8  0.00339924  …  8.01012e-9  4.78652e-9
 -9.36079e-8  0.0034061      8.3188e-9   4.24252e-9

[:, :, 4] =
  0.000307454  -0.0150898   -0.00391724  …   0.000301757
  6.05948e-5   -0.0152536   -0.00405883      8.45503e-5
  ⋮                                      ⋱
 -0.000307094  -0.00755439  -0.00156729  …  -0.000594673
 -0.000378125  -0.00746687  -0.00146262     -0.00051676
```

And the plot of these:

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.