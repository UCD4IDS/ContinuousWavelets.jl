# continuous 1-d; different scalings should lead to different sizes, different boundary condtions shouldn't
@testset "Construction Types" begin
    @testset "xSz=$xSize, b=$boundary, s=$s, d=$β, wfc=$(wfc.waveType)" for xSize = (33, 67), boundary = (DEFAULT_BOUNDARY, padded, NaivePer),
        s=[1, 2, 3.5, 8, 16], β = [1, 1.5,4.0],
        wfc in (wavelet(morl,s=s,boundary=boundary,
                        β=β),
                wavelet(dog0,s=s,boundary=boundary,β=β),
                wavelet(paul4,s=s,boundary=boundary,β=β),
                wavelet(cDb2, s=s, boundary=boundary,β=β))
        xc = rand(Float64,xSize)
        # the sizes are of course broken at this size, so no warnings needed
        yc = 3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            yc = cwt(xc, wfc)
        end
        if typeof(wfc.waveType) <: Union{Morlet, Paul}
            @test Array{ComplexF64, 2}==typeof(yc)
        else
            @test Array{Float64, 2}==typeof(yc)
        end
        nOctaves, totalWavelets, sRanges, sWidths =
            getNWavelets(xSize, wfc)
        @test totalWavelets == size(yc, 2)
    end
end
