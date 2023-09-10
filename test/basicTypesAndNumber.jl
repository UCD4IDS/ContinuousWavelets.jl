# continuous 1-d; different scalings should lead to different sizes, different boundary condtions shouldn't
@testset "Construction Types" begin
    xSizes = (33, 67, 128)
    boundaries = (DEFAULT_BOUNDARY, padded, NaivePer)
    sVals = [1, 2, 3.5, 8, 16]
    βs = [1, 1.5, 4.0]
    waveTypes = (morl, dog0, paul4, cDb2, morse)
    @testset "xSz=$xSize, b=$boundary, s=$s, β=$β, wfc=$(wave)" for xSize in xSizes,
        boundary in boundaries,
        s in sVals,
        β in βs,
        wave in waveTypes

        wfc = wavelet(wave, s = s, boundary = boundary, β = β)
        xr = rand(Float64, xSize)
        x32 = rand(Float32, xSize)
        xc = rand(ComplexF64, xSize)
        # the sizes are of course broken at this size, so no warnings needed
        yc = 3
        yr = 3
        y32 = 3
        with_logger(ConsoleLogger(stderr, Logging.Error)) do
            yr = cwt(xr, wfc)
            y32 = cwt(x32, wfc)
            yc = cwt(xc, wfc)
        end
        if isAnalytic(wfc.waveType)
            @test eltype(yr) == ComplexF64
            @test eltype(y32) == ComplexF32
            @test eltype(yc) == ComplexF64
        else
            @test eltype(yr) == Float64
            @test eltype(y32) == Float32
            @test eltype(yc) == ComplexF64
        end
        nOctaves, totalWavelets, sRanges, sWidths =
            getNWavelets(xSize, wfc)
        @test totalWavelets == size(yc, 2)
    end
end
# xSize = xSizes[3]; boundary = boundaries[3]; s = sVals[end]; β = βs[3]; wave = waveTypes[4]
# xSize = 33; boundary = SymBoundary(); s = 1; β = 1; wave = dog2
