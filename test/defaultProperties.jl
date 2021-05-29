wave = cwts[6]; bc = bcs[1]; β=βs[2]; ave = averagingLengths[1]; n = ns[1]
problemInds =Array{Any,1}()
@testset "Wavelet properties" begin
    bcs =(ContinuousWavelets.PerBoundary(),ContinuousWavelets.SymBoundary(), ContinuousWavelets.ZPBoundary())
    #cwts = (cDb2,ContinuousWavelets.dog1,ContinuousWavelets.dog2,ContinuousWavelets.paul1,ContinuousWavelets.morl,ContinuousWavelets.morse)
    cwts = (cDb2,ContinuousWavelets.Morse(3,2,1), ContinuousWavelets.Morse(1,1,1), ContinuousWavelets.Morse(1,4,1), ContinuousWavelets.Morse(2,3,1))
    βs =(1,4)
    ns = (128, 1382,2039)
    averagingLengths=(0,.5,1,2,4)
    @testset "length $n, with type $wave, bc $bc, β=$β, ave=$(ave)" for n in ns, wave in cwts, bc in bcs, β in βs, ave in averagingLengths
        wav = ContinuousWavelets.wavelet(wave, β=β,boundary=bc,averagingLength = ave)
        Ŵ,ω = (3,1)
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            Ŵ,ω = ContinuousWavelets.computeWavelets(n, wav)
        end
        # make sure the default last frequency isn't too large for the non-orthogonal derived wavelets
        if !(wave isa ContOrtho)
            t = @test minimum(abs.(Ŵ[end,:])/maximum(abs.(Ŵ)) .< 1e-3)
            if !(t isa Test.Pass)
                push!(problemInds,(wave,bc,β,ave,n))
            end
        end
        if wav.waveType isa Union{<:Morlet, <:Paul, <:Morse}
            nFreq, nSpace = ContinuousWavelets.setn(n, wav)
            fullŴ = [Ŵ; zeros(nFreq-1, size(Ŵ, 2))]
            fullŴ[end-nFreq+1:end,1] = reverse(Ŵ[:,1])
            spaceWaves = ifft(fullŴ,1)
        else
            if bc isa PerBoundary
                nLong = n
            elseif bc isa SymBoundary
                nLong = 2n
            elseif bc isa ZPBoundary
                nLong = 2^(ceil(Int,log2(n)))
            end
            spaceWaves = irfft(Ŵ,nLong,1)
        end
        # since the signals are centered at the boundary, the area that should be zero is in the middle
        nonSupported = spaceWaves[ceil(Int,n/2):end-ceil(Int,n/2),:]
        supported = spaceWaves[[1:floor(Int,n/2); end-floor(Int,n/2):end],:]
        # make sure the wavelets don't have significant support out into the mirrored signal
        if !(bc isa PerBoundary) && !(wave isa Paul) && size(Ŵ,2) > 1
            t = @test max([norm(nonSupported[:,i],Inf) / norm(supported[:,i],Inf) for i = 2:size(spaceWaves,2)]...) ≤ 1e-2
            if !(t isa Test.Pass)
                push!(problemInds,(wave,bc,β,ave,n))
            end
        elseif wave isa Paul && size(Ŵ,2) > 1 # the Paul wavelet decay is just too slow to avoid bleedover, so best to ignore the averaging size
            t = @test max([norm(nonSupported[:,i],Inf) / norm(supported[:,i],Inf) for i = 2:size(spaceWaves,2)]...) ≤ 1e-1
            if !(t isa Test.Pass)
                push!(problemInds,(wave,bc,β,ave,n))
            end
        end
        # if the averaging length is 0 and its only averaging at the sizes given, something is wrong
        t = @test size(Ŵ,2)>1 || ave>0
        if !(t isa Test.Pass)
            push!(problemInds,(wave,bc,β,ave,n))
        end
        if !(wave isa ContOrtho) && size(Ŵ,2) > 1
            # Guaranteed in a very different way for ContOrtho, and no guarantees made for just averaging
            t = @test max(abs.(Ŵ[end, :])...)  / norm(supported,Inf) ≤ 1e-1
            if !(t isa Test.Pass)
                push!(problemInds,(wave,bc,β,ave,n))
            end
        end
    end
end
# tools for manually inspecting what happened with the tests
if false
    wave,bc,β,ave,n = problemInds[3]
    wave = Morse(2,4,1)
    ave = 0
    wav = ContinuousWavelets.wavelet(wave, β=β,boundary=bc,averagingLength = ave)
    Ŵ,ω = ContinuousWavelets.computeWavelets(n, wav)
    plot(abs.(Ŵ), legend=false)
    ContinuousWavelets.getNthMoment(wav,1)^2
    ContinuousWavelets.getNthMoment(wav,2)
    ContinuousWavelets.normalizationConst(wav)
    plot(father(wav, ω, ContinuousWavelets.Father(),1))

    ContinuousWavelets.getMean(wav)
    ContinuousWavelets.getStd(wav)
    ContinuousWavelets.getNOctaves(128,wav)
    # what's the o.g. mother?
    nFreq, nSpace = ContinuousWavelets.setn(n, wav)
    Ŵ,ω = ContinuousWavelets.computeWavelets(128, wav)
    mo = mother(wav, 1,1,ω)
    fullMother = [mo; zeros(nFreq-1, size(mo, 2))]
    fullMother[end-nFreq+1:end,1] = reverse(mo[:,1])
    timeMom = ifft(fullMother,1)
    plot([real.(timeMom) imag.(timeMom)])
    if wav.waveType isa Union{<:Morlet, <:Paul, <:Morse}
        nFreq, nSpace = ContinuousWavelets.setn(n, wav)
        fullŴ = [Ŵ; zeros(nFreq-1, size(Ŵ, 2))]
        fullŴ[end-nFreq+1:end,1] = reverse(Ŵ[:,1])
        spaceWaves = ifft(fullŴ,1)
    else
        if bc isa PerBoundary
            nLong = n
        elseif bc isa SymBoundary
            nLong = 2n
        elseif bc isa ZPBoundary
            nLong = 2^(ceil(Int,log2(n)))
        end
        spaceWaves = irfft(Ŵ,nLong,1)
    end
    # since the signals are centered at the boundary, the area that should be zero is in the middle
    nonSupported = spaceWaves[ceil(Int,n/2):end-ceil(Int,n/2),:]
    supported = spaceWaves[[1:floor(Int,n/2); end-floor(Int,n/2):end],:]
    plot(real.(supported), legend=false)
    plot(real.(nonSupported), legend=false)
    plot(real.(spaceWaves), legend=false)
    getNWavelets(n,wav)
    ContinuousWavelets.getNOctaves(n,wav)
    ContinuousWavelets.getMinScaling(wav)
end
