```@meta ex
DocTestFilters = [
        r"\@ ContinuousWavelets .*",
        r"[ +-][0-9]\.[0-9]{3,5}e-1[5-9]",
        r"[ +-][0-9]\.[0-9]{3,5}e-[2-9][0-9]",
        r"im {2,7}",
    ]
```

# ContinuousWavelets

[![Build Status](https://travis-ci.com/dsweber2/ContinuousWavelets.jl.svg?branch=master)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
[![Coverage](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dsweber2.github.io/ContinuousWavelets.jl/dev/)

This package is an offshoot of [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) for the continuous wavelets.
Thanks to [Felix Gerick](https://github.com/fgerick) for the initial implementation there, with extension and further adaptation by David Weber and any other contributors listed on the right.
Currently, it implements 1D continuous wavelet transforms with the following mother wavelets:

![Mothers](https://dsweber2.github.io/ContinuousWavelets.jl/dev/mothers.svg)

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

julia> c = wavelet(Morlet(π), β=2);

julia> res = cwt(f, c)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.061863
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:7
2047×31 Matrix{ComplexF64}:
 -1.48637e-6+3.8241e-19im   …  0.000109978+9.67834e-5im
 -1.48602e-6+5.15534e-19im     -8.24922e-5+0.000130656im
            ⋮               ⋱             ⋮
 0.000435175+2.30636e-19im  …  -2.47195e-6-1.97048e-8im
 0.000435027-8.28725e-19im     -2.63499e-6+4.62331e-8im
```

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = cwt(f, c)

# output

2047×27 Matrix{Float64}:
 0.000926575  -0.00150445  …  -3.1189e-8   -2.84589e-8
 0.000926735  -0.00150491     -2.27135e-8  -1.96595e-8
 ⋮                         ⋱   ⋮
 0.000250708  -2.53146e-6  …  -3.11406e-9  -2.69593e-9
 0.000250698  -2.52694e-6     -4.29947e-9  -3.92694e-9
```

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = getMeanFreq(length(f), c)
dropped = icwt(res, c, DualFrames())
round.(dropped,sigdigits=12)

# output

┌ Warning: the canonical dual frame is off by 3.81e6, consider using one of the delta dual frames
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:41
2047-element Vector{Float64}:
 0.0069417253841
 0.00694223965946
 ⋮
 0.00268854944484
 0.00268852311536
```

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```jldoctest ex
julia> using Wavelets

julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0);

julia> res = circshift(cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.26753
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:12
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
  1.88768e-5   0.000266059  0.000196385  …   3.00171e-6
  8.31921e-5   0.000266939  0.000201693     -4.46452e-5
  ⋮                                      ⋱   ⋮            
  2.25286e-6   0.00198715   0.00125153   …   3.80696e-6
  2.64275e-6   0.0019801    0.00123431       3.47586e-6

[:, :, 2] =
  1.38846e-17  0.022676   0.00955726  …   4.66435e-18
  0.0          0.0226846  0.0095044      -2.60336e-18
  ⋮                                   ⋱   ⋮            
 -7.39019e-18  0.0341524  0.0108042   …  -1.86659e-18
  2.50852e-18  0.0342169  0.0107732      -1.48794e-18

[:, :, 3] =
 -4.25348e-7  0.00596895  0.00256794  …  1.86785e-8
 -4.37683e-7  0.00596787  0.0025637      7.83973e-9
  ⋮                                   ⋱  ⋮           
 -9.36725e-8  0.00339937  0.0016624   …  4.79906e-9
 -9.31697e-8  0.00340624  0.00165756     4.25503e-9

[:, :, 4] =
  0.000307892  -0.0150904   -0.0039183   …   0.000301771
  6.09684e-5   -0.0152542   -0.00405989      8.45675e-5
  ⋮                                      ⋱   ⋮            
 -0.000308175  -0.00755455  -0.00156652  …  -0.000594702
 -0.000378998  -0.00746703  -0.00146187     -0.000516786
```

And the plot of these:

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
