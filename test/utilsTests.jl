@testset "Coherence and Cross spectrum" begin
    rng = MersenneTwister(23425)
    Y = randn(rng, 2053, 4)
    c = wavelet(morl)
    YCohere = waveletCoherence(Y, Y, c)
    # test that the coherence of each signal with itself is approximately 1; note that if X is ever approximately zero in the scalogram this identity will break down for numerical reasons
    @test minimum([minimum(YCohere[:, :, ii, ii] .≈ 1) for ii = 1:size(YCohere, 3)])
    # test that every value is approximately at most 1
    @test minimum(YCohere .≤ 1 + sqrt(eps()))
    cYY = crossSpectrum(Y, Y, c)
    # cross spectrum of a signal with itself is an averaged wavelet power, so it should be real valued
    cYYDiag = cat([cYY[:, :, ii, ii] for ii = 1:size(cYY, 3)]..., dims = 3)
    @test minimum(real.(cYYDiag) .≈ abs.(cYYDiag))
end
