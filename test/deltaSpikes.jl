n = 2039;
n1 = n;
# wav = wavs[2]; k = ks[1]
@testset "Delta Spikes" begin
    wavs = (
        wavelet(cDb2),
        wavelet(cCoif4),
        wavelet(cBeyl),
    )
    ks = (1093, 408)
    @testset "at $k, with type $(wav.waveType)" for k in ks, wav in wavs
        x = zeros(n)
        x[k] = 1
        res = 3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            res = cwt(x, wav)
        end
        peaks = argmax(abs.(res), dims = 1)
        @test k - 2 <= peaks[end][1] <= k + 2 # give slight range to handle rounding issues
        Ŵ, ω = (3, 1)
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            Ŵ, ω = computeWavelets(n, wav)
        end
        W = irfft(Ŵ, 2n, 1)
        @test !any(abs.(circshift(W[:, 3], k - 1)[1:n] .- res[:, 3]) .> 1e-10)
        # make sure the wavelets don't have support out into the mirrored signal
    end
end
