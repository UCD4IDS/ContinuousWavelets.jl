using ContinuousWavelets, LinearAlgebra, FFTW,Wavelets, Interpolations
x = Wavelets.testfunction(n1, "Doppler"); plot(x)
W = wavelet(ContOrtho(wavelet(WT.coif2)), averagingLength=3, 
            β=2.0, normalization=2)
res = ContinuousWavelets.cwt(x, W)

# compare the convolution with a delta spike with continuous wavelet generated
# by Wavelets.jl



# all of the junk used to come up with ContOrtho in the first place
using MATLAB
mat"run('/home/dsweber/allHail/matlab/Wavelab850/WavePath.m')"
plot(mat"MakeWavelet(4,2^3,'Haar', 4, 'Mother',2048)'")
N = 256
mkWv = makeWavelet(wavelet(WT.db2),0); plot!(mkWv[3], [mkWv[1] mkWv[2]],legend=false); size(mkWv[2])
h =wavelet(WT.db2).qmf
using LinearAlgebra
[length(makinWavelet(wavelet(WT.db2),i)[1]) for i=0:12]
[length(makewavelet(wavelet(WT.db4),i)[1]) for i=0:12]
,labels= ["Mother" "Father"]
plot(abs.(rfft(mkWv[1])))
using Revise, Interpolations, Wavelets, ContinuousWavelets, FFTW
W = wavelet(ContOrtho(wavelet(WT.db2)),averagingLength=3,β=2.0)
c = W
typeof(W.waveType.o)
n1 = 1954
n,nSpace = ContinuousWavelets.setn(n1,c)
nOctaves, totalWavelets, sRange, sWidth = getNWavelets(n1, W)
isAve = (c.averagingLength > 0 && !(typeof(c.averagingType) <: NoAve)) ? 1 : 0
ω, daughters = ContinuousWavelets.analyticOrNot(c,n,totalWavelets)
d = ContinuousWavelets.findAveraging(c,ω,c.averagingType,φ,sRange,nSpace)
plot(d)
typeof(c)
daughters, ω = computeWavelets(n1,c)
plot(abs.(daughters),legend=false)
heatmap(abs.(daughters))

typeof(c.waveType )<:ContOrtho
N = 2*n1 # the total length
mS = 4 # minimum support size
w = wavelet(WT.db2)
qmf(W.waveType)
depth,ψLen = calcDepth(W,n1)
φ,ψ,t = makewavelet(W.waveType.o, depth) #\psi, \varphi
plot([φ ψ φ .* ψ])
sum(φ .* ψ)
itpφ = interpolate(ψ,BSpline(Quadratic(Reflect(OnGrid()))))
itpψ = interpolate(ψ,BSpline(Quadratic(Reflect(OnGrid()))))
using FFTW
abs.(rfft(φ))
abs.(rfft(ψ))
plot(ψ)
plot(abs.(rfft(φ)))
plot(abs.(rfft(ψ)))
N = n1
s = sRange[1]
subbedψ = itpψ(range(1,stop=ψLen,length =round(Int,N/s)))
subbedψ = itpψ(range(1,stop=ψLen,step = sRange[end]))
plot(subbedψ)
using Plots
plot(circshift(padTo(subbedψ, N),-round(Int,length(subbedψ)/2)))
plot(padTo(subbedψ, N))
daughters = zeros(N, length(sRange))
ψLen/N
N/16
for (ii,s) in enumerate(sRange)
    println("$N, $s, $(round(Int, N/s))")
    println("$(length(range(1, stop=ψLen, step= ψLen/N*s)))")
    daught = itpψ(range(1, stop=ψLen, step= ψLen/N *s))
    daughters[:, ii] = circshift(padTo(daught, N),
                                 -round(Int, length(daught)/2))
end
plot(ifftshift(daughters,1)[900:1050,:],legend=false)
plot(ifftshift(daughters,1)[:,:],legend=false)
log2.(sRange)
N*2
plot(daughters[:,7])
heatmap(abs.(rfft(daughters,1)))
plot(abs.(rfft(daughters[:,end-4:end], 1)))
plot(ifftshift(daughters[:,:],1))
daughts = rfft(daughters, 1);
daughts = daughts ./ [norm(x) for x in eachslice(daughts,dims=1)];
heatmap(abs.(daughts))
Morlet()
daughtersMorl, ω = computeWavelets(n1,wavelet(Morlet(), averagingLength=1))
wavelet(Morlet())
plot(ifftshift(irfft(daughtersMorl, n1, 1),1)[1:end, 2],legend=false)
1200-750
x = Wavelets.testfunction(n1, "Doppler"); plot(x)
W = wavelet(ContOrtho(wavelet(WT.coif2)), averagingLength=3, 
            β=2.0, normalization=2)
x = zeros(n1); x[round(Int,n1/2)] = 1
n1/2
res = ContinuousWavelets.cwt(x, W)
plot(res[:,1])
plot(res[:,2])
plot([x real.(res[:,1])])
heatmap(abs.(res))
plot(res[:,end-3])
typeof(wavelet(Morlet()))
nOctaves = log2(n1) - 4 - 2 + 1
using ContinuousWavelets
polySpacing(nOctaves, )
calcDepth(c, n1)
plot(Ψ)
heatmap(Ψ)
nOctaves, totalWavelets, sRange,sWidth = getNWavelets(N,wavelet(WT.morl))
sRange
nOctaves
2^nOctaves
sRange
Ψ
plot(abs.(daughts))
heatmap(abs.(daughts))
plot(daughters)
maximum(abs.(rfft(daughters,1))[:,end-1])
maximum(abs.(daughts)[:,end])
plot(abs.(rfft(daughters, 1)))
padTo(v, N) = cat(v, zeros(eltype(v), N-length(v)), dims = 1)
plot(subbedψ)
plot(ψ)
N/s
round(3.5)
wavelet(WT.db2).qmf
qmf(WT.db2)
"""
    nIters, sigLength = calcDepth(w,N)
given a wavelet type, calculate the depth necessary to get a
resulting wavelet longer than `N`. sigLength is the resulting length
"""
function calcDepth(w, N)
    origLen = length(qmf(w.waveType))
    netLen = origLen
    nIters = 0
    if N < netLen
        return nIters
    end
    while N > netLen - 2^nIters + 1
        netLen = nextLength(netLen, origLen)
        nIters += 1
    end
    return (nIters, netLen - 2^nIters + 1)
end
calcDepth(w,21)
nextLength(curr,orig) = curr*2 + orig - 1 # upsample then convolve
curr = length(w.qmf)
for i = 1:10
    global curr
    println(curr)
    curr = nextLength(curr,length(w.qmf))
end


makinWavelet(h, arg...) = makeWavelet(h.qmf, arg...)

function makinWavelet(h::AbstractVector, N::Integer=8)
    @assert N>=0
    sc = norm(h)
    h = h*sqrt(2)/sc
    phi = copy(h)
    psi = mirror(reverse(h))

    for i=1:N
        phi = DSP.conv(upsample(phi), h)
        psi = DSP.conv(upsample(psi), h)
    end
    return rmul!(phi,sc/sqrt(2)), rmul!(psi,sc/sqrt(2)), range(0, stop=length(h)-1, length=length(psi))
end
# comparing the wavelet computed by Wavelets.jl with the store ones made by computeWavelet.jl
# getting the right norm here is proving difficult, as they are somewhat different, which is to be expected
Ŵ, ω = computeWavelets(n, wav)
W = irfft(Ŵ, 2n,1)
ϕ,ψ,t = makewavelet(wav.waveType.o,6); size(ϕ)
plot(W[:,1])
plot(ifftshift(W[:,1])[1300:1700])
plot(ifftshift(W[abs.(W[:,1]) .> 1e-19,1]))
W[700,1]
n/174
2n/2^(wav.averagingLength)
abs.(W[30:end,1])
itpϕ = interpolate(ϕ,BSpline(Quadratic(Reflect(OnGrid()))))
firstLength = round(Int,2n/2^wav.averagingLength)
reϕ = itpϕ(range(1, stop = length(ϕ),length=firstLength))
justSupport = circshift(W[:,1], ceil(Int,firstLength/2))[1:firstLength]
(reϕ .- justSupport .* (norm(reϕ)/norm(justSupport)) ) ./ abs.(reϕ)
plot(reϕ ./ norm(reϕ) .- justSupport ./norm(justSupport))
plot([reϕ./norm(reϕ),justSupport ./norm(justSupport)])
plot(reϕ)
size(reϕ)
firstLength
plot(circshift(W[:,1], ceil(Int,firstLength/2)))
plot(circshift(W[:,1], firstLength)[1:firstLength])
ϕ,ψ,t = makewavelet(wav.waveType.o, 6)
size(ϕ)
wav = wavelet(cVaid, averagingLength = 2,normalization=2)
Ŵ,ω = computeWavelets(Int(length(ϕ)/2),wav)
W = irfft(Ŵ,length(ϕ),1)

plot(W,legend=false)
makewavelet()
cVaid
