# make a dataframe of the relative errors on a broad set of test cases
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
#wave = morl; bc = PerBoundary(); β=2; ave = -1; n = 1382; testF = "just Core"; inverseType = DualFrames()
# put the data frame in order of increasing error
sort!(df, :err)
# see the non-naive entries where the error is greater than 50%
show(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :], allrows=true)
sort!(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :], :dualErr)
show(last(df[(df.err .≥ .5) .& (map(x->x!=NaiveDelta(), df.invMeth)), :],50), allrows=true)
show(last(df[(df.N .== 128) .& (map(x->x==ZPBoundary(), df.bc)), :], 50), allrows=true)
show(last(df[(map(x->x==ZPBoundary(), df.bc)), :], 50), allrows=true)


inverseType=PenroseDelta()

βp = computeNonNegativeDualWeights(Ŵ)
heatmap(abs.(res)')
wav = wavelet(wave, β=β,boundary=bc,averagingLength = ave)
Ŵ,ω = computeWavelets(n,wav)
plot(Ŵ)
dual,dualErrr = ContinuousWavelets.getDualCoverage(n,wav,inverseType)
plot(dual)
Ŵdual = conj.(Ŵ ./ [norm(ŵ)^2 for ŵ in eachslice(Ŵ,dims=1)])
[norm(ŵ)^2 for ŵ in eachslice(Ŵ,dims=1)]
plot(Ŵ,legend=false)
plot(Ŵdual[:,1:end-0],legend=false)
res = cwt(x, wav)
xRecon = real.(ContinuousWavelets.icwt(res,wav,inverseType))
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
