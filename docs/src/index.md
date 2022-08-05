# ContinuousWavelets Documentation

Originally included in Wavelets.jl, this is a fork containing the types and methods specifically for doing continuous wavelet transforms. Current methods only include 1D wavelet transforms and their inverses.

The basic structure is similar to that of Wavelets.jl; first you choose one of the [Available Wavelet Families](@ref) of the `ContWaveClass` type, e.g. `Morlet(2π)`.
Then you set the general transform parameters via [CWT Construction](@ref), which specifies such properties as whether to average, the scaling rate, or the boundary conditions.
Finally, you perform the actual transform with`cwt`.

```@example basicEx
using Plots; gr(); #hide
Plots.reset_defaults(); #hide
using ContinuousWavelets, Plots, Wavelets, FFTW
n=2047;
f = testfunction(n, "Doppler");
p1=plot(f,legend=false,title="Doppler",xlims=(0,2000));
c = wavelet(Morlet(π), averagingType=NoAve(), β=2);
res = ContinuousWavelets.cwt(f, c)
p2=heatmap(abs.(res)', xlabel= "time index",
	ylabel="frequency index",colorbar=false);
l=@layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l);
```

There are 4 warnings you may see with some regularity:

```@example basicEx
julia> c = wavelet(Morlet(π), averagingType=NoAve(), β=2);

julia> res = ContinuousWavelets.cwt(f, c);
```

This comes from checking the value of the wavelets at 0 frequency.
For quasi-analytic wavelets such as Morlet wavelets, this means that there is significant non-zero mass in the negative frequency domain, which causes significant distortion.
Solvable either by increasing `aveLen` or increasing the mother wavelet mean frequency.

```@example basicEx
julia> c = wavelet(Morlet(1.5π), averagingType=NoAve(), β=2,extraOctaves=10);

julia> res = ContinuousWavelets.cwt(f, c);
```

This occurs because some of the constructed wavelets have significant mass beyond above the frequency resolution achievable for this signal length.
Usually solvable by simply decreasing `extraOctaves`.

```@example basicEx
julia> c = wavelet(Morlet(1.5π), averagingType=NoAve(), p=1);

julia> res = ContinuousWavelets.cwt(f, c);
```

These two almost always occur together as having peaks sufficiently far apart is an easy way for some frequencies to have insufficient coverage.
It is to some degree unavoidable for small values of `p` which governs which Fourier domain norm is preserved as we change the scale.

```@contents

```
