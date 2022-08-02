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
using ContinuousWavelets, Plots, Wavelets
n = 2047;
t = range(0,n/1000,length=n); # 1kHz sampling rate
f = testfunction(n, "Doppler");
p1 = plot(t, f,legend=false,title="Doppler",xticks=false)
c = wavelet(Morlet(π), β=2);
res = ContinuousWavelets.cwt(f, c)

# output

┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6
2047×31 Matrix{ComplexF64}:
 -2.84943e-6+3.86436e-19im    1.44519e-5-1.62298e-9im  …   0.000125192+0.000112017im   0.000109884+9.67013e-5im
 -2.84699e-6-6.11361e-20im    1.44509e-5-4.88515e-9im      -9.89185e-5+0.000149589im   -8.24222e-5+0.000130545im
 -2.84212e-6+4.37411e-20im    1.44489e-5-8.19594e-9im     -0.000180219-7.17018e-5im   -0.000153333-5.64666e-5im
 -2.83483e-6-3.11387e-19im    1.44459e-5-1.15878e-8im       2.75325e-5-0.000205094im    1.90839e-5-0.00016841im
 -2.82514e-6-1.31096e-19im    1.44419e-5-1.5093e-8im       0.000217137-3.06784e-5im    0.000172696-2.56466e-5im
 -2.81306e-6-3.38731e-19im    1.44369e-5-1.87442e-8im  …   0.000101677+0.000208577im    7.79501e-5+0.000162848im
 -2.79865e-6-9.8753e-19im     1.44309e-5-2.25737e-8im     -0.000166245+0.0001764im    -0.000128919+0.000132755im
 -2.78192e-6+4.91292e-20im    1.44239e-5-2.66139e-8im     -0.000231546-8.56946e-5im   -0.000172323-6.71036e-5im
 -2.76293e-6+5.80924e-19im    1.44159e-5-3.08974e-8im      -1.72186e-5-0.000244897im   -9.39619e-6-0.000179998im
            ⋮                                          ⋱                                          ⋮
 0.000172941+2.7571e-19im   -0.000580229+4.73334e-5im      -2.59937e-6-7.00815e-7im    -2.59966e-6-7.0039e-7im
 0.000171274+1.41585e-19im  -0.000580627+4.17758e-5im      -2.58041e-6-6.22483e-7im    -2.58362e-6-6.24181e-7im
 0.000169814-7.90531e-21im  -0.000580975+3.62141e-5im  …   -2.56582e-6-5.32728e-7im    -2.56497e-6-5.41108e-7im
 0.000168561-5.81895e-20im  -0.000581273+3.06488e-5im      -2.57532e-6-4.46165e-7im    -2.56135e-6-4.4822e-7im
 0.000167516-3.07438e-19im  -0.000581522+2.50805e-5im        -2.584e-6-4.03348e-7im    -2.57801e-6-3.83539e-7im
 0.000166679-6.64104e-19im  -0.000581721+1.95096e-5im      -2.51699e-6-3.66859e-7im    -2.54192e-6-3.54776e-7im
 0.000166051-1.45091e-18im  -0.000581871+1.39368e-5im       -2.4225e-6-2.1293e-7im      -2.4431e-6-2.37279e-7im
 0.000165633+5.67188e-19im   -0.00058197+8.36266e-6im  …   -2.48306e-6+6.63823e-9im    -2.46986e-6-1.96881e-8im
 0.000165423+1.25225e-18im   -0.00058202+2.78765e-6im      -2.65444e-6+5.18306e-8im    -2.63276e-6+4.61939e-8im

```

And now we make a scalogram to actually interpret these entries:

```jldoctest ex
freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
freqs[1] = 0
p2 = heatmap(t, freqs, log.(abs.(res).^2)', xlabel= "time (s)", ylabel="frequency (Hz)", colorbar=false, c=cgrad(:viridis, scale=:log10))
l = @layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)

# output
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6
Plot{Plots.GRBackend() n=2}
```

![Doppler](/docs/doppler.svg)

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```jldoctest ex
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
c = wavelet(dog2, β = 2);
res = ContinuousWavelets.cwt(f, c)

# output

```

```jldoctest ex
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = ContinuousWavelets.icwt(res, c, DualFrames())

# output

```

```jldoctest ex
p1 = plot(f, legend = false, title = "Smoothing and dropping bumps", linewidth = 2)
plot!(dropped, linewidth = 3)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout = l)
# output
```

![Bumps](/docs/bumps.svg)

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```julia
exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims = 2)
c = wavelet(cDb2, β = 2, extraOctaves = -0);
res = circshift(ContinuousWavelets.cwt(exs, c), (0,1,0))
```

![parallel transforms](/docs/multiEx.svg)

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
- Various additional wavelet families, such as Morse wavelets.
