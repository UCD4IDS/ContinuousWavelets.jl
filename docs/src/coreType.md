# Wavelet types
There are two tiers of wavelet types in this package. The most abstract is the
`ContWave` type, representing a class of wavelets
```@docs
ContWave
```

```@setup basicEx
using ContinuousWavelets, Plots, Wavelets, FFTW, Logging
using Plots; gr()
Plots.reset_defaults()
global_logger(Logging.SimpleLogger(stderr,Logging.Error))
n= 2047;
function mapTo(waveType, isReal=true,window=2047-100:2047+100;d=2)
	if isReal
		c = wavelet(waveType,β=d)
		morlet,ω = computeWavelets(n,c)
		return ifftshift(irfft(morlet,2*n,1))[window,:]
	else
		c = wavelet(waveType,β=d)
		morlet,ω = computeWavelets(n,c)
		return ifftshift(ifft(morlet,1))[window,:]
	end
end
tmp = mapTo(Morlet(π),false, 1024-50:1024+50)[:,2]
p1=plot([real.(tmp) imag.(tmp)],title="Morlet",
	labels=["real" "imaginary"],ticks=nothing,linewidth=5);
tmp = mapTo(paul2,false, 1024-50:1024+50)[:,4]
p2=plot([real.(tmp) imag.(tmp)],title="Paul 2",
	labels=["real" "imaginary"],ticks=nothing,linewidth=5);
p3=plot(mapTo(dog2)[:,2],title="derivative of gaussians (dog2)",legend=false,ticks=nothing,linewidth=5);
p4=plot(mapTo(cHaar)[:,2],title="Haar",legend=false,ticks=nothing,linewidth=5);
p5=plot(mapTo(cBeyl;d=1)[:,2],title="Beylkyin",legend=false,ticks=nothing,linewidth=5);
p6=plot(mapTo(cVaid;d=1)[:,2],title="Vaidyanthan",legend=false,ticks=nothing,linewidth=5);
p7=plot(mapTo(cDb2;d=1)[:,2],title="Daubhechies 2",legend=false,ticks=nothing,linewidth=5);
p8=plot(mapTo(cCoif2;d=1)[:,2],title="Coiflet 2",legend=false,ticks=nothing,linewidth=5);
p9=plot(mapTo(cSym4;d=1)[:,2],title="Symlet 4",legend=false,ticks=nothing,linewidth=5);
p10=plot(mapTo(cBatt4;d=1)[:,2],title="Battle-Lemarie, 4",legend=false,ticks=nothing,linewidth=5);
plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,layout=(2,5),size=300 .*(5, 2.2));
savefig("mothers.svg")#hide
```
![](mothers.svg)
Above are examples of every mother wavelet defined in this package; the only
analytic and/or complex wavelets are the `Morlet` and the `Paul` wavelets. Once
you have chosen a type of wavelet, this is used to construct the more specific
CWT, which specifies more details of the transform, such as frequency spacing,
whether to include an averaging filter or not, a frame upper bound, etc.
