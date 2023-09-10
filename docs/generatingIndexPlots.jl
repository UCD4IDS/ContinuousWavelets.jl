include("make.jl")
doctest(ContinuousWavelets)
using ContinuousWavelets, Plots, Wavelets
n = 2047;
t = range(0, n / 1000, length = n); # 1kHz sampling rate
f = testfunction(n, "Doppler");
p1 = plot(t, f, legend = false, title = "Doppler", xticks = false)
c = wavelet(Morlet(π), β = 2);
res = ContinuousWavelets.cwt(f, c)
freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
freqs[1] = 0
p2 = heatmap(
    t,
    freqs,
    log.(abs.(res) .^ 2)',
    xlabel = "time (s)",
    ylabel = "frequency (Hz)",
    colorbar = false,
    c = cgrad(:viridis, scale = :log10),
)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout = l)
savefig("doppler.svg")



f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2);
res = ContinuousWavelets.cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
p2 = heatmap(
    1:n,
    freqs,
    abs.(res)',
    xlabel = "time (ms)",
    ylabel = "Frequency (Hz)",
    colorbar = false,
    c = :viridis,
)
dropped = ContinuousWavelets.icwt(res, c, DualFrames())
p1 = plot(f, legend = false, title = "Smoothing and dropping bumps", linewidth = 2)
plot!(dropped, linewidth = 3)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout = l)
savefig("bumps.svg")



exs = cat(
    testfunction(n, "Doppler"),
    testfunction(n, "Blocks"),
    testfunction(n, "Bumps"),
    testfunction(n, "HeaviSine"),
    dims = 2,
);
c = wavelet(cDb2, β = 2, extraOctaves = -0)
res = circshift(ContinuousWavelets.cwt(exs, c), (0, 1, 0))

p1 = plot(
    plot(
        exs[:, 1],
        legend = false,
        title = "Doppler",
        yticks = [],
        xticks = [],
        linewidth = 2,
    ),
    plot(
        exs[:, 2],
        legend = false,
        title = "Blocks",
        yticks = [],
        xticks = [],
        linewidth = 2,
    ),
    plot(
        exs[:, 3],
        legend = false,
        title = "Bumps",
        yticks = [],
        xticks = [],
        linewidth = 2,
    ),
    plot(
        exs[:, 4],
        legend = false,
        title = "HeaviSine",
        yticks = [],
        xticks = [],
        linewidth = 2,
    ),
    layout = (1, 4),
)
p2 = plot(
    heatmap(
        identity.(res[:, :, 1])',
        xticks = false,
        yticks = [],
        c = :viridis,
        colorbar = false,
    ),
    heatmap(
        identity.(res[:, :, 2])',
        xticks = false,
        yticks = [],
        c = :viridis,
        colorbar = false,
    ),
    heatmap(
        identity.(res[:, :, 3])',
        xticks = false,
        yticks = [],
        c = :viridis,
        colorbar = false,
    ),
    heatmap(
        identity.(res[:, :, 4])',
        xticks = false,
        yticks = [],
        c = :viridis,
        colorbar = false,
    ),
    layout = (1, 4),
)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout = l)
savefig("multiEx.svg")
