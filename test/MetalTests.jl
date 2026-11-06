using BenchmarkTools

@testset "GPU vs CPU (Metal)" begin
    metal_available = @isdefined(Metal) && Metal.functional()
    if metal_available
        @testset "reflect - Metal" begin
            boundaries = (PerBoundary(), SymBoundary(), ZPBoundary())
            ns = (37, 200, 256)  # 200 is not a power of 2; 256 already is (padding = 0)
            @testset "n1=$n1, boundary=$bt" for n1 in ns, bt in boundaries
                x_cpu = randn(Float32, n1, 3)  # extra trailing dim: reflect must handle N > 1
                x_mtl = MtlArray(x_cpu)

                r_cpu = reflect(x_cpu, bt)
                r_mtl = reflect(x_mtl, bt)

                @test r_mtl isa MtlArray
                @test eltype(r_mtl) == eltype(x_mtl)
                @test size(r_mtl) == size(r_cpu)
                @test Array(r_mtl) ≈ r_cpu
            end
        end

        waveTypes = (morl, dog2, paul2, Morse(3, 20, 1))
        β = 2
        boundaries = (PerBoundary(), SymBoundary(), ZPBoundary())
        averagingLength = 2
        extraOctaves = 0
        xSizes = (2048, 8192, 32768)

        @testset "xSz=$xSize, boundary=$boundary, wave=$wave" for xSize in xSizes,
            boundary in boundaries,
            wave in waveTypes

            wfc = wavelet(wave, β = β, boundary = boundary,
                averagingLength = averagingLength,
                extraOctaves = extraOctaves)

            # CPU arrays:
            xr_cpu = randn(Float32, xSize)
            xc_cpu = randn(ComplexF32, xSize)

            # Metal arrays:
            xr_mtl = MtlArray(xr_cpu)
            xc_mtl = MtlArray(xc_cpu)

            # Compute daughters on CPU, then move to Metal:
            daughters_cpu, ω = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                computeWavelets(xSize, wfc)
            end
            daughters_mtl = MtlArray(Float32.(daughters_cpu))

            # CPU transforms:
            yr_cpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xr_cpu, wfc, Float32.(daughters_cpu))
            end
            yc_cpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xc_cpu, wfc, copy(Float32.(daughters_cpu)))
            end

            # Metal transforms:
            yr_mtl = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xr_mtl, wfc, daughters_mtl)
            end
            yc_mtl = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xc_mtl, wfc, MtlArray(Float32.(daughters_cpu)))
            end

            # Correctness: Check that the Metal result matches the CPU result.
            @test Array(yr_mtl) ≈ convert.(eltype(yr_mtl), yr_cpu)
            @test Array(yc_mtl) ≈ convert.(eltype(yc_mtl), yc_cpu)

            # Make sure that element types are preserved.
            @test eltype(Array(yr_mtl)) <: Union{Float32,ComplexF32}
            @test eltype(Array(yc_mtl)) <: Union{Float32,ComplexF32}
            @test eltype(yr_cpu) <: Union{Float32,ComplexF32}
            @test eltype(yc_cpu) <: Union{Float32,ComplexF32}

            # Compare the speed:
            daughters_cpu_bench = Float32.(daughters_cpu)
            t_cpu = @belapsed cwt($xr_cpu, $wfc, copy($daughters_cpu_bench))
            t_mtl = @belapsed begin
                cwt($xr_mtl, $wfc, copy($daughters_mtl))
                Metal.synchronize()
            end

            @info "CWT speed: wave=$(wave), xSize=$(xSize), boundary=$(boundary)" cpu_time=t_cpu mtl_time=t_mtl speedup=t_cpu/t_mtl
        end

        @testset "device mismatch handling" begin
            xSize = 512
            wfc = wavelet(morl, β = 2)
            x_cpu = randn(Float32, xSize)
            x_mtl = MtlArray(x_cpu)

            daughters_cpu, _ = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                computeWavelets(xSize, wfc)
            end
            daughters_cpu = Float32.(daughters_cpu)
            daughters_mtl = MtlArray(daughters_cpu)

            @test_throws "same device" cwt(x_mtl, wfc, daughters_cpu)

            y_mtl = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(x_mtl, wfc, daughters_mtl)
            end
            @test Array(icwt(y_mtl, wfc, PenroseDelta())) ≈
                  icwt(Array(y_mtl), wfc, PenroseDelta())
            @test Array(icwt(y_mtl, wfc, NaiveDelta())) ≈
                  icwt(Array(y_mtl), wfc, NaiveDelta())
            @test Array(icwt(y_mtl, wfc, DualFrames())) ≈
                  icwt(Array(y_mtl), wfc, DualFrames())
        end

        @testset "Float64 rejection" begin
            wfc = wavelet(morl, β = 2)
            @test_throws Exception cwt(MtlArray(randn(Float64, 256)), wfc)
        end
    else
        @warn "Metal not available or not functional, skipping GPU tests"
    end
end