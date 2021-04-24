#wave = cwts[4]; bc = bcs[1]; β=βs[2]; ave = averagingLengths[end]; n = ns[1]; testF = typesOfTestFunctions[1]; inverseType = inversionMethods[3]
#wave = morl; bc = ZPBoundary(); ave = -1; n = 1382; testF = "just Core"; inverseType = DualFrames(); β = 1.5; eOct = 0
@testset "Inversion" begin
    bcs =(ContinuousWavelets.PerBoundary(), ContinuousWavelets.ZPBoundary(),)
    cwts = (ContinuousWavelets.dog2,ContinuousWavelets.morl)
    βs =(1,)
    averagingLengths=(0,)
    extraOctaves=(0,)
    typesOfTestFunctions = ["Doppler"]
    inversionMethods = [NaiveDelta(), PenroseDelta(), DualFrames()]
    ns = (128, 2039)
    @testset "length $n, with type $wave, bc $bc, β=$β, ave=$(ave) ex=$(testF), inv=$(inverseType)" for n in ns, wave in cwts, bc in bcs, β in βs, ave in averagingLengths, testF in typesOfTestFunctions, inverseType in inversionMethods, eOct in extraOctaves
        if !(wave isa Morlet && eOct < 0)
            if testF =="just Core"
                x̂ = zeros(ComplexF64, n>>1+1); x̂[n>>2 .+ (-10:10)] = 5*(11.5 .- abs.(range(-10/2,10/2,length=21))) .* randn(ComplexF64, 21)
                x = irfft(x̂, n) # test function that is definitely solidly in the "safe" range
            else
                x = testfunction(n,testF)
            end

            wav = ContinuousWavelets.wavelet(wave, β=β,boundary=bc,averagingLength = ave, extraOctaves = eOct)
            with_logger(ConsoleLogger(stderr, Logging.Error)) do
                res = ContinuousWavelets.cwt(x, wav)
                xRecon = real.(ContinuousWavelets.icwt(res, wav, inverseType))
                err = norm(xRecon-x)/norm(x)
                @test err < 3 # none of them are wildly off
            end
        end
    end
end
