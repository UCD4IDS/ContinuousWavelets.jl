using ContinuousWavelets, Plots, Wavelets, FFTW,Logging
global_logger(SimpleLogger(min_level=Logging.Error))
global_logger(Logging.SimpleLogger(stderr,Logging.Error))
n=2047;
f = testfunction(n, "Doppler");
p1=plot(f,legend=false,title="Doppler",xlims=(0,2000))
c = wavelet(Morlet(π), averagingType=NoAve(), β=6);
res = cwt(f, c)
p2=heatmap(abs.(res)', xlabel= "time index", ylabel="frequency index",colorbar=false)
l=@layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l);
savefig("doppler.svg");#hide
DEFAULT_BOUNDARY; supertype(typeof(morl))
a,b=computeWavelets(n,c)
heatmap(a)
using LaTeXStrings
dRate = 4
waveType = Morlet()
Ψ1 = wavelet(waveType, s=8, β =dRate, averagingLength=2)
pyplot()
locs = ContinuousWavelets.polySpacing(8,Ψ1);
scatter(1:length(locs), locs, legend=:bottomright, label="mean log frequency",
        xlabel="Wavelet Index (x)", ylabel= "log-Frequency", color=:black)
scatter!(length(locs):length(locs),locs[end:end],markersize=10,markershape=:x,color=:black, label=:none)
firstW, lastW,stepS = ContinuousWavelets.genSamplePoints(8,Ψ1)
b = (dRate/Ψ1.Q)^(1 ./dRate)*(8+Ψ1.averagingLength)^((dRate-1)/dRate)
t= range(1,stop=length(locs),step=.1)
curve = b .*(range(firstW,stop=(locs[end]/b)^dRate,length=length(t))).^(1 / dRate)
plot!(t, curve, color=:blue, line=:dash, label=L"b(x+\gamma)^{^1/_\beta}",
      legend=:bottomright,
      xrange=(0,length(locs)+3), xticks= [1; 5:5:1+length(locs)...],
      yrange=(minimum(locs)-1, maximum(locs)+1),
      yticks=(2:2:10,["Ave.Length", (4:2:8)..., "N.Octaves"]))
x = range(15, stop=28, step=.5)
y(x)= curve[end] .+ b/dRate*24 .^(1/dRate-1).*(x .-24)
plot!(x, y(x), c=:black,line=2,label=:none)
annotate!(length(locs)-.125, locs[end]+7/16, 
          Plots.text(L"\frac{dy}{dx}=^{1}/_{Q}", 9, :black, :center))
savefig("plotOfLogCentralFrequencies.pdf")

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
     colorbar_title=L"\widehat{\psi_i}")#hide

savefig("changeBeta.svg")
aasdf = randn(10,10)
plot(aasdf, aspect_ratio=1.0)
plot(aasdf)
l=@layout [grid(4,3) [a{.5w};b{.5w};c{.5w};d{.5w}]]
ps = [heatmap(i*j*rand(10,10),colorbar=false) for i=1:4,j=1:4]
plot(ps...,layout=l,size=(1100,800))
Q=8
eval(:Q)
