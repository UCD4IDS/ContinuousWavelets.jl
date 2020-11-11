# Wavelet Frequency spacing #
Frequently, using a fixed ratio for scaling the wavelets results in too many
large scale wavelets. There are several ways of dealing with this; in this
package, the scaling factors have the form $2^{a(mx+x_0)^{^1/_\beta}}$, for
suitable choice of $a$,$m$, $x_0$, and $\beta$.  The figure gives an example
of the chosen scaling factors in log frequency. 

```@setup waves
using ContinuousWavelets, Plots, Wavelets, FFTW, LaTeXStrings, Logging
global_logger(Logging.SimpleLogger(stderr,Logging.Error))
dRate = 4
waveType = Morlet()
Ψ1 = wavelet(waveType, s=8, β =dRate, averagingLength=2)
pyplot()
locs = ContinuousWavelets.polySpacing(8,Ψ1);
scatter(1:length(locs), locs, legend=:bottomright, label="mean log frequency",
        xlabel="Wavelet Index (x)", ylabel= "log-Frequency (y)", color=:black)
scatter!(length(locs):length(locs),locs[end:end],markersize=10,markershape=:x,color=:black, label=:none)
firstW, lastW,stepS = ContinuousWavelets.genSamplePoints(8,Ψ1)
b = (dRate/Ψ1.Q)^(1 ./dRate)*(8+Ψ1.averagingLength)^((dRate-1)/dRate)
t= range(1,stop=length(locs),step=.1)
curve = b .*(range(firstW,stop=(locs[end]/b)^dRate,length=length(t))).^(1 / dRate)
plot!(t, curve, color=:blue, line=:dash, label=L"y=a(mx+x_0)^{^1/_\beta}",
      legend=:bottomright,
      xrange=(0,length(locs)+3), xticks= [1; 5:5:1+length(locs)...],
      yrange=(minimum(locs)-1, maximum(locs)+1),
      yticks=(2:2:10,["Ave.Length", (4:2:8)..., "N.Octaves"]));
x = range(15, stop=28, step=.5)
y(x)= curve[end] .+ b/dRate*24 .^(1/dRate-1).*(x .-24)
plot!(x, y(x), c=:black,line=2,label=:none);
annotate!(length(locs)-.125, locs[end]+7/16, 
          Plots.text(L"\frac{dy}{dx}=^{1}/_{Q}", 9, :black, :center));
savefig("plotOfLogCentralFrequencies.svg")
```
![](plotOfLogCentralFrequencies.svg)

If $\beta$ is 1, then we have a linear relation between the index and the
log-frequency, and $Q$ gives exactly the number of wavelets per octave
throughout. As $\beta$ increases, the wavelets skew more and more heavily to
high frequencies. The default value is 4.

The user chooses $\beta$, $Q$ (the number of wavelets per octave at the last
point), and Ave. Length (the number of octaves covered by the averaging
function), and then $a$, $m$, $x_0$, and the total number of wavelets $N_w$ are
chosen so that:
1. The first wavelet is scaled by $2^{\textrm{Ave. Length}}$, so the curve
    $a(mx+x_0)^{^1/_\beta}$ goes through the point 
    $(x,y)=(1,\textrm{Ave. Length})$.
2. The derivative $\frac{\mathrm{d}y}{\mathrm{d}x}$ at the last point is
    $\frac{1}{Q}$, so the "instantaneous" number of wavelets $x$ per octave $y$
    is $Q$. Each type of wavelet has a maximum scaling $2^{N_{Octaves}}$
    returned by `getNOctaves` (generally half the signal length), so the final
    point $N_w$ satisfies both $y(N_w) = N_{Octaves}$ and
    $y'(N_w)=^1/_Q$.
3. Finally, the spacing is chosen so that there are exactly $Q$ wavelets in the
   last octave.

If you are interested in the exact computation, see the function `polySpacing`.
As some examples of how the wavelet bank changes as we change $\beta$:
```@example waves
n=2047
Ψ1 = wavelet(morl, s=8, β=1)
d1, ξ = computeWavelets(n,Ψ1)
Ψ2 = wavelet(morl, s=8, β =2)
d2, ξ = Wavelets.computeWavelets(n,Ψ2)
Ψ4 = wavelet(morl, s=8, β =4)
d4, ξ = Wavelets.computeWavelets(n,Ψ4)
matchingLimits = (minimum([d1 d2 d4]), maximum([d1 d2 d4]))# for 
plot(heatmap(1:size(d1,2), ξ, d1, color=:Greys, 
             yaxis = (L"\omega", ), 
             xaxis = ("wavelet index", ),
             title=L"\beta=1"*" ("*L"\Psi1"*")", colorbar=false,
             clims=matchingLimits),
     heatmap(1:size(d2,2), ξ, d2, color=:Greys, 
             yticks=[], 
             xaxis = ("wavelet index", ),
             title=L"\beta=2"*" ("*L"\Psi2"*")", colorbar=false,
             clims=matchingLimits), 
     heatmap(1:size(d4,2), ξ, d4,color=:Greys, yticks=[],
             xaxis = ("wavelet index", ),
             title=L"\beta=4"*" ("*L"\Psi4"*")"),
     layout=(1,3),  clims=matchingLimits, 
     colorbar_title=L"\widehat{\psi_i}");#hide
savefig("changeBeta.png") #hide
```
![](changeBeta.png)

note that the low-frequency coverage increases drastically as we decrease
$\beta$.
