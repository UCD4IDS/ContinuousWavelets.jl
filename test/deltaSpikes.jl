n = 2039
@testset "Delta Spikes" begin
    @testset "at $k, with type $(wav.waveType)" for k in (1093, 408), wav in (wavelet(cDb2), wavelet(cCoif4), wavelet(cBeyl))
        x = zeros(n); x[k]=1
        res=3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            res = cwt(x,wav)
        end
        peaks = argmax(abs.(res), dims=1)
        @test peaks[end][1] == k
        Ŵ,ω = (3,1)
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            Ŵ,ω = computeWavelets(n,wav)
        end
        W = irfft(Ŵ, 2n,1)
        @test !any(abs.(circshift(W[:,1],k-1)[1:n] .- res[:,1]) .> 1e-10)
    end
end
