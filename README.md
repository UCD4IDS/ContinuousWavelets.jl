# ContinuousWavelets

[![Build Status](https://github.com/UCD4IDS/ContinuousWavelets.jl/actions/workflows/CI.yml/badge.svg)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
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

```julia
julia> using ContinuousWavelets, Plots, Wavelets

julia> n = 2047;

julia> t = range(0, n / 1000, length=n); # 1kHz sampling rate

julia> f = testfunction(n, "Doppler");

julia> p1 = plot(t, f, legend=false, title="Doppler", xticks=false)
Plot{Plots.PyPlotBackend() n=1}

julia> c = wavelet(Morlet(π), β=2)
CWT{Morlet mean 3.141592653589793, Father Wavelet, Q=8.0, β=2.0,aveLen=0.0, frame=1.0, norm=Inf, extraOctaves=0.0}

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

And now we make a scalogram to actually visualize these entries:

```julia
freqs = getMeanFreq(computeWavelets(n, c)[1])
freqs[1] = 0
p2 = heatmap(t, freqs, log.(abs.(res).^2)', xlabel= "time (s)", ylabel="frequency (Hz)", colorbar=false, c=cgrad(:viridis, scale=:log10))
l = @layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)
```

![Doppler](/docs/doppler.svg)

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```julia
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
c = wavelet(dog2, β = 2)
res = cwt(f, c);
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = getMeanFreq(length(f), c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = icwt(res, c, DualFrames());
p1 = plot(f, legend=false, title="Smoothing and dropping bumps", linewidth=2)
plot!(dropped, linewidth=3)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)
```

![Bumps](/docs/bumps.svg)

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```julia
julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0)
CWT{Continuous db2, Father Wavelet, Q=8.0, β=2.0,aveLen=0.0, frame=1.0, norm=Inf, extraOctaves=0.0}

julia> res = circshift(cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.26753
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:10
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
 1.88768e-5  0.000266059  …  4.67422e-5   3.00171e-6
 8.31921e-5  0.000266939     1.56546e-5  -4.46452e-5
 ⋮                        ⋱  ⋮
 2.25286e-6  0.00198715   …  4.24046e-6   3.80696e-6
 2.64275e-6  0.0019801       4.37922e-6   3.47586e-6

[:, :, 2] =
  1.82235e-17  0.022676   0.00955726  …   3.2542e-18
 -1.73557e-18  0.0226846  0.0095044      -4.77282e-18
  ⋮                                   ⋱
  1.53669e-18  0.0341524  0.0108042   …   6.05043e-19
  7.9713e-18   0.0342169  0.0107732       7.14045e-19

[:, :, 3] =
 -4.25348e-7  0.00596895  …  4.4713e-8   1.86785e-8
 -4.37683e-7  0.00596787     3.30062e-8  7.83973e-9
  ⋮                       ⋱  ⋮
 -9.36725e-8  0.00339937  …  7.99443e-9  4.79906e-9
 -9.31697e-8  0.00340624     8.3032e-9   4.25503e-9

[:, :, 4] =
  0.000307892  -0.0150904   -0.0039183   …   0.000301771
  6.09684e-5   -0.0152542   -0.00405989      8.45675e-5
  ⋮                                      ⋱
 -0.000308175  -0.00755455  -0.00156652  …  -0.000594702
 -0.000378998  -0.00746703  -0.00146187     -0.000516786
```

And the plot of these:

```julia
p1 = plot(plot(exs[:, 1], legend=false, title="Doppler", yticks=[], xticks=[], linewidth=2), plot(exs[:, 2], legend=false, title="Blocks", yticks=[], xticks=[], linewidth=2), plot(exs[:, 3], legend=false, title="Bumps", yticks=[], xticks=[], linewidth=2), plot(exs[:, 4], legend=false, title="HeaviSine", yticks=[], xticks=[], linewidth=2), layout=(1, 4))
p2 = plot(heatmap(identity.(res[:, :, 1])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 2])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 3])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 4])', xticks=false, yticks=[], c=:viridis, colorbar=false), layout=(1, 4))
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)
```

![parallel transforms](/docs/multiEx.svg)

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
