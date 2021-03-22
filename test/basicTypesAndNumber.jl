# continuous 1-d; different scalings should lead to different sizes, different boundary condtions shouldn't
@testset "Construction Types" begin
    xSizes = (33, 67, 128)
    boundaries = (DEFAULT_BOUNDARY, padded, NaivePer)
    sVals = [1, 2, 3.5, 8, 16]
    βs = [1, 1.5,4.0]
    waveTypes =(morl, dog0, paul4, cDb2)
    @testset "xSz=$xSize, b=$boundary, s=$s, β=$β, wfc=$(wave)" for xSize in xSizes, boundary in boundaries, s in sVals, β in βs, wave in waveTypes
        wfc = wavelet(wave, s=s, boundary=boundary,β=β)
        xc = rand(Float64, xSize);
        # the sizes are of course broken at this size, so no warnings needed
        yc = 3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            yc = cwt(xc, wfc);
        end
        if typeof(wfc.waveType) <: Union{Morlet, Paul}
            @test Array{ComplexF64, 2}==typeof(yc)
        else
            @test Array{Float64, 2}==typeof(yc)
        end
        nOctaves, totalWavelets, sRanges, sWidths = getNWavelets(xSize, wfc);
        @test totalWavelets == size(yc, 2)
    end
end
# xSize = xSizes[1]; boundary = boundaries[1]; s = sVals[1]; β = βs[1]; wave = waveTypes[1]
