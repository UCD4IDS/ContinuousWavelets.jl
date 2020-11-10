# ContinuousWavelets Documentation #
Originally included in Wavelets.jl, this is a fork containing the types and
methods specifically for doing continuous wavelet transforms. Current methods
only include 1D wavelet transforms and their inverses.

The basic structure is similar to that of Wavelets.jl; first you choose one of
the ![Wavelet types](@ref) of the `ContWaveClass` type, e.g. `Morlet(2π)`.
Then you set the general transform parameters by constructing a ![CWT
type](@ref), which specifies such properties as whether to average,
the scaling rate, or the boundary conditions. Finally, you perform the actual
transform with`cwt`.

```@example basicEx
using ContinuousWavelets, Plots, Wavelets, FFTW
using Logging#hide
global_logger(Logging.SimpleLogger(stderr,Logging.Error))#hide
n=2047;
f = testfunction(n, "Doppler");
p1=plot(f,legend=false,title="Doppler",xlims=(0,2000));
c = wavelet(Morlet(π), averagingType=NoAve(), β=2);
res = cwt(f, c)
p2=heatmap(abs.(res)', xlabel= "time index", 
	ylabel="frequency index",colorbar=false)
l=@layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)
savefig("doppler.svg");#hide
```
![](doppler.svg)

```@contents
```
