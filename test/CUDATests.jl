using BenchmarkTools

@testset "GPU vs CPU (CUDA)" begin
    cuda_available = @isdefined(CUDA) && CUDA.functional()
    if cuda_available
        @testset "reflect - CUDA" begin
            boundaries = (PerBoundary(), SymBoundary(), ZPBoundary())
            ns = (37, 200, 256)  # 200 is not a power of 2; 256 already is (padding = 0)
            @testset "n1=$n1, boundary=$bt" for n1 in ns, bt in boundaries
                x_cpu = randn(Float32, n1, 3)  # extra trailing dim: reflect must handle N > 1
                x_gpu = CuArray(x_cpu)

                r_cpu = reflect(x_cpu, bt)
                r_gpu = reflect(x_gpu, bt)

                @test r_gpu isa CuArray
                @test eltype(r_gpu) == eltype(x_gpu)
                @test size(r_gpu) == size(r_cpu)
                @test Array(r_gpu) ≈ r_cpu
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

            # GPU arrays:
            xr_gpu = CuArray(xr_cpu)
            xc_gpu = CuArray(xc_cpu)

            # Compute daughters on CPU, then move to GPU:
            daughters_cpu, ω = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                computeWavelets(xSize, wfc)
            end
            daughters_gpu = CuArray(Float32.(daughters_cpu))

            # CPU transforms:
            yr_cpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xr_cpu, wfc, Float32.(daughters_cpu))
            end
            yc_cpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xc_cpu, wfc, copy(Float32.(daughters_cpu)))
            end

            # GPU transforms:
            yr_gpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xr_gpu, wfc, daughters_gpu)
            end
            yc_gpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(xc_gpu, wfc, CuArray(Float32.(daughters_cpu)))
            end

            # Correctness: Check that the GPU result matches the CPU result.
            @test Array(yr_gpu) ≈ convert.(eltype(yr_gpu), yr_cpu)
            @test Array(yc_gpu) ≈ convert.(eltype(yc_gpu), yc_cpu)

            # Make sure that element types are preserved.
            @test eltype(Array(yr_gpu)) <: Union{Float32,ComplexF32}
            @test eltype(Array(yc_gpu)) <: Union{Float32,ComplexF32}
            @test eltype(yr_cpu) <: Union{Float32,ComplexF32}
            @test eltype(yc_cpu) <: Union{Float32,ComplexF32}

            # Compare the speed:
            daughters_cpu_bench = Float32.(daughters_cpu)
            t_cpu = @belapsed cwt($xr_cpu, $wfc, copy($daughters_cpu_bench))
            t_gpu = @belapsed begin
                cwt($xr_gpu, $wfc, copy($daughters_gpu))
                CUDA.synchronize()
            end

            @info "CWT speed: wave=$(wave), xSize=$(xSize), boundary=$(boundary)" cpu_time=t_cpu gpu_time=t_gpu speedup=t_cpu/t_gpu
        end

        @testset "device mismatch handling" begin
            xSize = 512
            wfc = wavelet(morl, β = 2)
            x_cpu = randn(Float32, xSize)
            x_gpu = CuArray(x_cpu)

            daughters_cpu, _ = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                computeWavelets(xSize, wfc)
            end
            daughters_cpu = Float32.(daughters_cpu)
            daughters_gpu = CuArray(daughters_cpu)

            # cwt with a GPU signal but CPU daughters (an easy mistake: forgot
            # to move `daughters` over) should fail fast with a clear message,
            # not the opaque GPU-kernel-compilation KernelError this used to
            # produce.
            @test_throws "same device" cwt(x_gpu, wfc, daughters_cpu)

            # icwt on a fully-GPU forward transform should now work
            # transparently for all three InverseType variants: β/
            # canonDualFrames are built in res's precision and adapted onto
            # res's device before being combined with it.
            y_gpu = with_logger(ConsoleLogger(stderr, Logging.Error)) do
                cwt(x_gpu, wfc, daughters_gpu)
            end
            @test Array(icwt(y_gpu, wfc, PenroseDelta())) ≈ icwt(Array(y_gpu), wfc, PenroseDelta())
            @test Array(icwt(y_gpu, wfc, NaiveDelta()))   ≈ icwt(Array(y_gpu), wfc, NaiveDelta())
            @test Array(icwt(y_gpu, wfc, DualFrames()))   ≈ icwt(Array(y_gpu), wfc, DualFrames())
        end
    else
        @warn "CUDA not available or not functional, skipping GPU tests"
    end
end