```@setup invEx
using ContinuousWavelets, Plots, Wavelets, FFTW, Logging
using Plots; gr()
Plots.reset_defaults()
global_logger(Logging.SimpleLogger(stderr,Logging.Error))
```
# Wavelet Inversion #
The continuous wavelet transform is a redundant shift-invariant frame transform. As such, there isn't a single inverse transform, although there is a canonical pseudo-inverse. For more see, for example, chapter 5 of [A Wavelet Tour of Signal Processing](https://wavelet-tour.github.io/).

In this package, we have 3 pseudo-inverses:
```@docs
icwt
```
