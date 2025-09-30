# Available Wavelet Families
There are two tiers of wavelet types in this package. The most abstract is the `ContWave` type, representing a class of wavelets.
This is split into several strictly continuous wavelets, and a `ContOrtho<:ContWave` type, which is a supertype of continuous versions of the orthogonal wavelets defined in [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl).
```@docs
ContWave
```

```@setup basicEx
using ContinuousWavelets, Plots, Wavelets, FFTW, Logging
using Plots; gr()
Plots.reset_defaults()
n= 2047;
function mapTo(waveType, isReal = true, window = 1:2047; d = 1, γ = 4.0, β = 2.0, cf = 1.0, kwargs...)
    if waveType == Morse
        morse_wav = Morse(float(γ), float(β), float(cf))
        c = wavelet(morse_wav; kwargs...)
    else
        c = wavelet(waveType; β=d, kwargs...)
    end
    waves, ω = computeWavelets(n, c)
    if isReal
        return circshift(irfft(waves, 2*n, 1), (1024, 0))[window, :]
    else
        waves = cat(waves, zeros(2047, size(waves, 2)), dims = 1)
        return circshift(ifft(waves, 1), (1024, 0))[window, :]
    end
end
tmp = mapTo(Morlet(π), false; averagingLength = -0.2)[:, 2]
p1 = plot([real.(tmp) imag.(tmp)], title = "Morlet", labels = ["real" "imaginary"], ticks = nothing, linewidth = 5)
tmp = mapTo(paul2, false, averagingLength = -0.5)[:, 2]
p2 = plot([real.(tmp) imag.(tmp)], title = "Paul 2", labels = ["real" "imaginary"], ticks = nothing, linewidth = 5)
tmpMorse1 = mapTo(Morse, false; β=3, γ=10.0, cf=1.0, averagingLength=-1)[:, 2]
p3 = plot([real.(tmpMorse1) imag.(tmpMorse1)], title = "Morse (β=3, γ=10)", labels = ["real" "imaginary"], ticks = nothing, linewidth = 4)
tmpMorse2 = mapTo(Morse, false; β=1, γ=3.0, cf=1.0, averagingLength=-2)[:, 2] 
p4 = plot([real.(tmpMorse2) imag.(tmpMorse2)], title = "Morse (β=1, γ=3)", labels = ["real" "imaginary"], ticks = nothing, linewidth = 4)
p5 = plot(mapTo(dog2; averagingLength = -1.5)[:, 2], title = "derivative of gaussians (dog2)", legend = false, ticks = nothing, linewidth = 5)
p6 = plot(mapTo(cHaar, true; averagingLength = 1)[:, 2], title = "Haar", legend = false, ticks = nothing, linewidth = 5)
p7 = plot(mapTo(cBeyl, true; d = 1, averagingLength = -0)[:, 2], title = "Beylkyin", legend = false, ticks = nothing, linewidth = 5)
p8 = plot(mapTo(cVaid, true; d = 1, averagingLength = -0)[:, 2], title = "Vaidyanathan", legend = false, ticks = nothing, linewidth = 5)
p9 = plot(mapTo(cDb2; d = 1, averagingLength = -0)[:, 2], title = "Daubhechies 2", legend = false, ticks = nothing, linewidth = 5)
p10 = plot(mapTo(cCoif2, true; d = 1, averagingLength = -0)[:, 2], title = "Coiflet 2", legend = false, ticks = nothing, linewidth = 5)
p11 = plot(mapTo(cSym4, true; d = 1, averagingLength = -0)[:, 2], title = "Symlet 4", legend = false, ticks = nothing, linewidth = 5)
k = 0600;
p12 = plot(mapTo(cBatt4, true, 1024-k:1024+k; d = 1, averagingLength = -1)[:, 2], title = "Battle-Lemarie, 4", legend = false, ticks = nothing, linewidth = 5);
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout = (3, 4), size = 300 .* (5, 3.2))
savefig("mothers.svg")#hide
```
![](mothers.svg)
Above are examples of every mother wavelet family defined in this package; the only analytic and/or complex wavelets are the `Morlet` and the `Paul` wavelet families.
Once you have chosen a type of wavelet, this is used to construct the more specific CWT, which specifies more details of the transform, such as frequency spacing, whether to include an averaging filter or not, a frame upper bound, etc.
