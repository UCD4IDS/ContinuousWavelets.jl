x = testfunction(n,"Doppler") # this is maybe a uniquely difficult example due to high freq stuff
x = testfunction(n,"Bumps")
125.5 .- abs.(range(-251/2,251/2,length=251))
x̂ = zeros(ComplexF64, 1020); x̂[350:600] = 5*(125.5 .- abs.(range(-251/2,251/2,length=251))) .* randn(ComplexF64, 251)
x̂ = zeros(ComplexF64, 1020); x̂[250] = 1
x = irfft(x̂, n) # test function that is definitely solidly in the "safe" range
plot(x)
plot(abs.(rfft(x)))
plot(x)
sum(x)
@views plot(sum(pinv(Ŵ[1:lastReasonableFreq,:]),dims=2))
@views lastReasonableFreq = findlast(abs.(Ŵ[:,end])/maximum(abs.(Ŵ[:,end])).>.5)
plot(Ŵ[1:1478, :], labels=false)
size(Ŵ)
2^(round(Int, log2(1382)))
zero(0)
wave = cwts[4]; bc = bcs[1]; β=βs[2]; ave = averagingLengths[end]; n = ns[1]; testF = typesOfTestFunctions[1]; inverseType = inversionMethods[3]
wave = morl; bc = ZPBoundary(); β=2; ave = -1; n = 1382; testF = "just Core"; inverseType = PenroseDelta()
using DataFrames
df = DataFrame(err = Float64[], dualErr=Float64[], fnc=String[], bc=WaveletBoundary[], cwt = ContinuousWavelets.ContWaveClass[], β = Int[], aveLen = Int[], invMeth = ContinuousWavelets.InverseType[], N = Int[])
@testset "Inversion" begin
    bcs =(PerBoundary(),SymBoundary(), ZPBoundary())
    cwts = (cDb2,dog1,dog2,paul1,morl)
    βs =(1, 2)
    averagingLengths=(-1,0,1)
    typesOfTestFunctions = ["Doppler" "Bumps" "Blocks" "HeaviSine" "just Core"]
    inversionMethods = [NaiveDelta(), PenroseDelta(), DualFrames()]
    ns = (128, 1382,2039)
    @testset "length $n, with type $wave, bc $bc, β=$β, ave=$(ave) ex=$(testF)" for n in ns, wave in cwts, bc in bcs, β in βs, ave in averagingLengths, testF in typesOfTestFunctions, inverseType in inversionMethods
        if testF =="just Core"
            x̂ = zeros(ComplexF64, n>>1+1); x̂[n>>2 .+ (-10:10)] = 5*(11.5 .- abs.(range(-10/2,10/2,length=21))) .* randn(ComplexF64, 21)
            x = irfft(x̂, n) # test function that is definitely solidly in the "safe" range
        else
            x = testfunction(n,testF)
        end
        if inverseType isa NaiveDelta
            errTol = 3 # naiveDelta really isn't very good at all
        elseif inverseType isa DualFrames && wave isa Morlet
            errTol = 1 # generally good unless the frames are particularly bad at the end
        else
            errTol = .01 # everything else should have error less than 1%
        end

        wav = wavelet(wave, β=β,boundary=bc,averagingLength = ave)
        res = 3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            res = cwt(x, wav)
            xRecon = real.(ContinuousWavelets.icwt(res, wav, inverseType))
            err = norm(xRecon-x)/norm(x)
            dualWav,dualErr = ContinuousWavelets.getDualCoverage(n, wav, DualFrames());
            push!(df, (err, dualErr, testF, bc, wave, β, ave, inverseType, n))
            #@test  err < errTol
        end
    end
end
# put the data frame in order of increasing error
sort!(df, :err)
# see the non-naive entries where the error is greater than 50%
show(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :], allrows=true)
sort!(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :], :dualErr)
show(last(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :],50), allrows=true)
x = zeros(n>>1+1); x[100:250].=1; x = irfft(x,n)
res = cwt(x,wav)
xRecon = ContinuousWavelets.icwt(res,wav)
plot([x real.(xRecon)])
plot(abs.(rfft(real.(xRecon) + irfft(im*rfft(imag.(xRecon)),size(xRecon,1)))))
plot(plot([real.(rfft(real.(xRecon))) imag.(rfft(real.(xRecon)))],labels=["real" "imag"],title="xRecon"), plot([real.(rfft(x)) imag.(rfft(x))],labels=["real" "imag"],title="x"), plot([real.(rfft(real.(xRecon)-x)) imag.(rfft(real.(xRecon)-x))],labels=["real" "imag"],title="Diff"))
plot(xRecon)
norm(x-xRecon)/norm(x)
dual, dualNorm = ContinuousWavelets.getDualCoverage(n, wav)
plot(abs.(dual))
res = cwt(x, wav)
plot(sum(abs.(ContinuousWavelets.explicitConstruction(Ŵ)).^2,dims=2).^(1/2))
plot(abs.(ContinuousWavelets.explicitConstruction(Ŵ)),legend=false)
Ŵ, ω = computeWavelets(n,wav)
plot(abs.(Ŵ),legend=false)
heatmap(abs.(res)')
β = ContinuousWavelets.computeDualWeights(Ŵ, wav); dualWav,dualCover = ContinuousWavelets.getDualCoverage(n, wav, PenroseDelta());
plot(plot([real.(β') imag.(β')],labels=["real" "imag"],title="β"), plot([real.(dualWav) imag.(dualWav)],title="dual total"),labels=["real" "imag"])
dualCover
w = 20; plot([real.(Ŵ[:,w]) imag.(Ŵ[:,w])],title="dual total")
dualCover
plot(real.(dualWav))
tmp = [3+im; zeros(199)]
plot([real.(fft(tmp)) imag.(fft(tmp))])
plot(irfft(im*ones(100),199))
plot(real.(ifft([0 im*ones(100); -im*ones(100)])))
plot(real.(xRecon))
plot([abs.(rfft(x)) abs.(rfft(xRecon))])
plot(abs.(rfft(x)))
        println()
using NonNegLeastSquares
n1= n
plot(sum(β .* Ŵ,dims=2))
plot(dual)
plot(sum(Ŵ .* netMultiply, dims=2))
plot(Ŵ .* netMultiply)
plot(sum(Ŵ .* netMultiply,dims=2))
xRecon = real.(sum(res .* β, dims=2))
plot(abs.(rfft(xRecon)))
plot!(abs.(rfft(x)))
function computeNonNegativeDualWeights(Ŵ)
    @views lastReasonableFreq = argmax(abs.(Ŵ[:,end]))
    β = nonneg_lsq(Ŵ[1:lastReasonableFreq,:], ones(lastReasonableFreq))
    Wdag = pinv(Ŵ[1:lastReasonableFreq,:])
    #@views Wdag[:,2:end] .*=2
    β = sum(Wdag, dims=2)
    β[2:end] .*= 2
    return β
end
# the paul wavelets are, as expected, severely ill conditioned b/c of the extremely slow decay
βp = computeNonNegativeDualWeights(Ŵ)
wav = wavelet(morl;β=2,p=p,boundary=PerBoundary())
heatmap(abs.(res)')
Ŵ,ω = computeWavelets(n,wav)
plot(Ŵ)
dual = ContinuousWavelets.getDualCoverage(n,wav)
plot(dual)
res = cwt(x, wav)
xRecon = real.(ContinuousWavelets.icwt(res,wav))
plot([x xRecon .*norm(x)/norm(xRecon)],labels=["original" "reconstructed"])
plot([x xRecon],labels=["original" "reconstructed"])
norm(xRecon), norm(x)
norm(x-xRecon .*norm(x)/norm(xRecon))
norm(x-xRecon)/norm(x)
plot(imag(xRecon))
plot(abs.(xRecon))
plot(abs.(fft(xRecon)))
plot([x xRecon],labels=["original" "reconstructed"])
plot([x x-xRecon .*norm(x)/norm(xRecon)],labels=["original" "error"])
plot([x x-xRecon],labels=["original" "error"])
plot([abs.(rfft(x)) abs.(rfft(xRecon))],labels=["original" "reconstruction"])
plot([abs.(rfft([x; reverse(x)])) abs.(rfft([xRecon; reverse(xRecon)]))],labels=["original" "reconstruction"])
plot([abs.(rfft(x)) abs.(rfft(x-xRecon .*norm(x)/norm(xRecon)))],labels=["original" "error"])
plot(plot(x), plot(xRecon))
plot(plot(abs.(rfft(x))), plot(abs.(rfft(xRecon))))
plot(plot(abs.(rfft([x; reverse(x)]))), plot(abs.(rfft(xRecon))))
plot([x xRecon])
plot([x x-xRecon .*norm(x)/norm(xRecon)])
