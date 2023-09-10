# Wavelet Frequency Spacing

Frequently, using a fixed ratio for scaling the wavelets results in too many large scale wavelets.
There are several ways of dealing with this; in this package, the scaling factors have the form $2^{x_0 +mx^{^1/_\beta}}$, for suitable choice of $a$,$m$, $x_0$, and $\beta$.
The user chooses $\beta$ and $Q$, and the rest are chosen to match the derivative at the last frequency to be $^{1}/_{Q}$, as in the figure.

```@setup waves
using ContinuousWavelets, Plots, Wavelets, FFTW, LaTeXStrings, Logging
using Plots;
gr();
Plots.reset_defaults();
dRate = 4
waveType = Morlet()
Q = 8
nOct = 8
Ψ1 = wavelet(waveType, s=Q, β=dRate, averagingLength=2)
# sketch of how the frequencies are chosen
locs = polySpacing(nOct, Ψ1)
a = getMinScaling(Ψ1) + Ψ1.averagingLength
β = Ψ1.β
lastWavelet = Ψ1.Q * (nOct - a)
b = 1 / Q * lastWavelet^((β - 1) / β)
t = range(1, stop=length(locs), step=0.1)
curve = a .+ b .* (range(0, stop=lastWavelet, length=length(t))) .^ (1 / β)
x = range(12, stop=22, step=0.5)
ycord(x) = locs[end] .+ b / β * lastWavelet .^ (1 / β - 1) .* (x .- length(locs)) * (lastWavelet - 1) / length(locs)
y = ycord.(x)
#Figure 3.1
scatter(1:length(locs), locs, legend=:bottomright, label="mean log frequency", xlabel="Wavelet Index (x)", ylabel="log-Frequency (y)", color=:black)
scatter!(length(locs):length(locs), locs[end:end], markersize=8, markershape=:x, color=:black, facecolor=:black, label=:none)
plot!(t, curve, color=:blue, line=:dash, label=L"y=x_0 + mx^{^1/_\beta}", legend=:bottomright, legendfontsize=12, xrange=(0, length(locs) + 3), xticks=[1; 5:5:1+length(locs)...], yrange=(minimum(locs) - 1, maximum(locs) + 1), yticks=(range(floor(Int, minimum(locs)), ceil(Int, maximum(locs)), step=2), [L"\alpha", (4:2:6)..., "N.Octaves"]))
plot!(x, y, color=:black, line=1.5, label=:none)
annotate!(length(locs) - 1 / 8, locs[end] + 7.5 / 16, Plots.text(L"\frac{\mathrm{d}y}{\mathrm{d}x}=^{1}/_{Q}", 11, :black, :center))
savefig("plotOfLogCentralFrequencies.svg")
```

![](plotOfLogCentralFrequencies.svg)

If $\beta$ is 1, then we have a linear relation between the index and the log-frequency, and $Q$ gives exactly the number of wavelets per octave throughout.
As $\beta$ increases, the wavelets skew more and more heavily to high frequencies.
The default value is 4.

The user chooses `β` (the frequency decay), `Q` (the number of wavelets per octave at the last point), and `aveLen` (the number of octaves covered by the averaging function, here $x_0$), and then $m$, the number of wavelets $N_w$, and the spacing of $x$ are chosen so that:

1. The derivative $\frac{\mathrm{d}y}{\mathrm{d}x}$ at the last point is $^{1}/_{Q}$, so the "instantaneous" number of wavelets $x$ per octave $y$ is $Q$.
   Each type of wavelet has a maximum scaling $2^{N_{Octaves}}$ returned by `getNOctaves` (generally half the signal length), so the final point $N_w$ satisfies both $y(N_w) = N_{Octaves}$ and $y'(N_w)=^1/_Q$.
2. The spacing is chosen so that there are exactly $Q$ wavelets in the last octave.

If you are interested in the exact computation, see the function `polySpacing`.
As some examples of how the wavelet bank changes as we change $\beta$:

```@example waves

n=2047
Ψ1 = wavelet(morl, s=8, β=1)
d1, ξ = computeWavelets(n,Ψ1)
Ψ2 = wavelet(morl, s=8, β =2)
d2, ξ = computeWavelets(n,Ψ2)
Ψ4 = wavelet(morl, s=8, β =4)
d4, ξ = computeWavelets(n,Ψ4)
matchingLimits = (minimum([d1 d2 d4]), maximum([d1 d2 d4]))#hide
plot(heatmap(1:size(d1,2), ξ, d1, color=:Greys, yaxis = (L"\omega", ), xaxis = ("wavelet index", ), title=L"\beta=1"*" ("*L"\Psi1"*")", colorbar=false, clims=matchingLimits),  heatmap(1:size(d2,2), ξ, d2, color=:Greys, yticks=[], xaxis = ("wavelet index", ), title=L"\beta=2"*" ("*L"\Psi2"*")", colorbar=false, clims=matchingLimits),  heatmap(1:size(d4,2), ξ, d4,color=:Greys, yticks=[], xticks=[1, 5, 10, 14, 18], xaxis = ("wavelet index", ), title=L"\beta=4"*" ("*L"\Psi4"*")"), layout=(1,3), clims=matchingLimits, colorbar_title=L"\widehat{\psi_i}")
savefig("changeBeta.png") #hide
```

![](changeBeta.png)

note that the low-frequency coverage increases drastically as we decrease
$\beta$.
