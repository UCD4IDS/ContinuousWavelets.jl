# ContinuousWavelets

[![Build Status](https://travis-ci.com/dsweber2/ContinuousWavelets.jl.svg?branch=master)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
[![Coverage](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dsweber2.github.io/ContinuousWavelets.jl/dev/)

This package is an offshoot of [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) for the continuous wavelets.
Thanks to [Felix Gerick](https://github.com/fgerick) for the initial implementation there, with extension and further adaptation by David Weber and any other contributors listed on the right.
Currently, it implements 1D continuous wavelet transforms with the following mother wavelets:

![Mothers](/docs/MotherWavelets.svg)

Which covers several standard continuous wavelet families, both real and analytic, as well as continuous versions of the orthogonal wavelet transforms implemented in [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl).

Basic Usage
---------
Install via the package manager and load with `using`

```julia
julia> Pkg.add("Wavelets")
julia> using Wavelets
```

Basic usage example on a doppler test function. 
```julia
using ContinuousWavelets, Plots, Wavelets
n=2047;
t = range(0,n/1000,length=n); # 1kHz sampling rate
f = testfunction(n, "Doppler");
p1=plot(t, f,legend=false,title="Doppler",xticks=false)
c = wavelet(Morlet(π), β=2);
res = cwt(f, c)
# plotting
freqs = getMeanFreq(computeWavelets(n, c)[1])
freqs[1] = 0
p2=heatmap(t,freqs, abs.(res)', xlabel= "time (s)", ylabel="frequency (Hz)",colorbar=false)
l=@layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)
```
![Doppler](/docs/doppler.svg)

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDual()`, `PenroseDual()`, and `DualFrame()`. As a toy example:

``` julia
f = testfunction(n, "Bumps");
p1=plot(f,legend=false,title="Bumps",xlims=(0,2000))
c = wavelet(dog2, β=2);
res = cwt(f, c)
# dropping the middle peaks
res[620:1100,:] .=0
# and smoothing the remaining peaks
res[:,10:29] .= 0
# actually doing the inversion
dropped = ContinuousWavelets.icwt(res,c,DualFrames())

p1 = plot([dropped f],legend=false, title="Smoothing and dropping bumps")
p2=heatmap(abs.(res)', xlabel= "time index", ylabel="frequency index",colorbar=false)
l=@layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)
```
![Bumps](/docs/bumps.svg)

It can also handle collections of examples at the same time, should you need to do a batch of transforms:
``` julia
exs = cat(testfunction(n, "Doppler"),testfunction(n,"Blocks"),testfunction(n,"Bumps"),testfunction(n,"HeaviSine"),dims=2)
c = wavelet(cDb2, β=2,extraOctaves=-0);
res = cwt(exs, c)
```
![parallel transforms](/docs/multiEx.svg)
There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

Possible extensions
------------
- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
- Various additional wavelet families, such as Morse wavelets.
