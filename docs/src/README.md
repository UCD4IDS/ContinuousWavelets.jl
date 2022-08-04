# ContinuousWavelets

[![Build Status](https://travis-ci.com/dsweber2/ContinuousWavelets.jl.svg?branch=master)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
[![Coverage](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dsweber2.github.io/ContinuousWavelets.jl/dev/)

This package is an offshoot of [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) for the continuous wavelets.
Thanks to [Felix Gerick](https://github.com/fgerick) for the initial implementation there, with extension and further adaptation by David Weber and any other contributors listed on the right.
Currently, it implements 1D continuous wavelet transforms with the following mother wavelets:

![Mothers](docs/mothers.svg)

Which covers several standard continuous wavelet families, both real and analytic, as well as continuous versions of the orthogonal wavelet transforms implemented in [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl).

## Basic Usage

Install via the package manager and load with `using`

```julia
julia> Pkg.add("ContinuousWavelets")
julia> using ContinuousWavelets
```

Basic usage example on a doppler test function.

```jldoctest ex
julia> using ContinuousWavelets, Plots, Wavelets

julia> n = 2047;

julia> t = range(0, n / 1000, length=n); # 1kHz sampling rate

julia> f = testfunction(n, "Doppler");

julia> p1 = plot(t, f, legend=false, title="Doppler", xticks=false)
Plot{Plots.PyPlotBackend() n=1}

julia> c = wavelet(Morlet(π), β=2)

julia> res = ContinuousWavelets.cwt(f, c)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6
2047×31 Matrix{ComplexF64}:
 -2.84943e-6+3.86436e-19im    1.44519e-5-1.62298e-9im  …   0.000125192+0.000112017im   0.000109884+9.67013e-5im
 -2.84699e-6-6.11361e-20im    1.44509e-5-4.88515e-9im      -9.89185e-5+0.000149589im   -8.24222e-5+0.000130545im
 -2.84212e-6+4.37411e-20im    1.44489e-5-8.19594e-9im     -0.000180219-7.17018e-5im   -0.000153333-5.64666e-5im
 -2.83483e-6-3.11387e-19im    1.44459e-5-1.15878e-8im       2.75325e-5-0.000205094im    1.90839e-5-0.00016841im
 -2.82514e-6-1.31096e-19im    1.44419e-5-1.5093e-8im       0.000217137-3.06784e-5im    0.000172696-2.56466e-5im
 -2.81306e-6-3.38731e-19im    1.44369e-5-1.87442e-8im  …   0.000101677+0.000208577im    7.79501e-5+0.000162848im
 -2.79865e-6-9.8753e-19im     1.44309e-5-2.25737e-8im     -0.000166245+0.0001764im    -0.000128919+0.000132755im
 -2.78192e-6+4.91292e-20im    1.44239e-5-2.66139e-8im     -0.000231546-8.56946e-5im   -0.000172323-6.71036e-5im
 -2.76293e-6+5.80924e-19im    1.44159e-5-3.08974e-8im      -1.72186e-5-0.000244897im   -9.39619e-6-0.000179998im
            ⋮                                          ⋱                                          ⋮
 0.000172941+2.7571e-19im   -0.000580229+4.73334e-5im      -2.59937e-6-7.00815e-7im    -2.59966e-6-7.0039e-7im
 0.000171274+1.41585e-19im  -0.000580627+4.17758e-5im      -2.58041e-6-6.22483e-7im    -2.58362e-6-6.24181e-7im
 0.000169814-7.90531e-21im  -0.000580975+3.62141e-5im  …   -2.56582e-6-5.32728e-7im    -2.56497e-6-5.41108e-7im
 0.000168561-5.81895e-20im  -0.000581273+3.06488e-5im      -2.57532e-6-4.46165e-7im    -2.56135e-6-4.4822e-7im
 0.000167516-3.07438e-19im  -0.000581522+2.50805e-5im        -2.584e-6-4.03348e-7im    -2.57801e-6-3.83539e-7im
 0.000166679-6.64104e-19im  -0.000581721+1.95096e-5im      -2.51699e-6-3.66859e-7im    -2.54192e-6-3.54776e-7im
 0.000166051-1.45091e-18im  -0.000581871+1.39368e-5im       -2.4225e-6-2.1293e-7im      -2.4431e-6-2.37279e-7im
 0.000165633+5.67188e-19im   -0.00058197+8.36266e-6im  …   -2.48306e-6+6.63823e-9im    -2.46986e-6-1.96881e-8im
 0.000165423+1.25225e-18im   -0.00058202+2.78765e-6im      -2.65444e-6+5.18306e-8im    -2.63276e-6+4.61939e-8im


```

And now we make a scalogram to actually interpret these entries:

```jldoctest ex; output=false
freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
freqs[1] = 0
p2 = heatmap(t, freqs, log.(abs.(res).^2)', xlabel= "time (s)", ylabel="frequency (Hz)", colorbar=false, c=cgrad(:viridis, scale=:log10))
l = @layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)

# output
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:6
Plot{Plots.PyPlotBackend() n=2}
```

![Doppler](/docs/doppler.svg)

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```jldoctest ex; output = false
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c)

# output

2047×27 Matrix{Float64}:
 0.00262067   -0.00150482  -2.16134e-6  -9.32656e-7  -6.06247e-7  …  -3.75559e-8   -3.42124e-8   -3.11967e-8   -2.84659e-8
 0.00262075   -0.00150528  -2.16503e-6  -9.32531e-7  -6.05646e-7     -2.985e-8     -2.61034e-8   -2.27191e-8   -1.96643e-8
 0.00262089   -0.00150621  -2.17246e-6  -9.32285e-7  -6.0445e-7      -1.90126e-8   -1.53614e-8   -1.22211e-8   -9.55694e-9
 0.00262111   -0.0015076   -2.18368e-6  -9.31925e-7  -6.02672e-7     -9.98394e-9   -7.26844e-9   -5.16474e-9   -3.5861e-9
 0.00262141   -0.00150945  -2.1988e-6   -9.31463e-7  -6.00331e-7     -4.71452e-9   -3.17363e-9   -2.13149e-9   -1.45404e-9
 0.00262177   -0.00151176  -2.21794e-6  -9.30914e-7  -5.97451e-7  …  -2.45418e-9   -1.71933e-9   -1.26651e-9   -9.85913e-10
 0.00262221   -0.00151453  -2.24129e-6  -9.30297e-7  -5.94064e-7     -1.73806e-9   -1.36213e-9   -1.11251e-9   -9.3166e-10
 0.00262272   -0.00151775  -2.26904e-6  -9.29635e-7  -5.90205e-7     -1.58516e-9   -1.31798e-9   -1.11305e-9   -9.48609e-10
 0.00262331   -0.00152143  -2.30144e-6  -9.28953e-7  -5.85916e-7     -1.58454e-9   -1.33761e-9   -1.13556e-9   -9.66579e-10
 ⋮                                                                ⋱                               ⋮
 0.000849844  -2.6919e-6   -1.4791e-7   -9.09068e-8  -6.38853e-8     -1.16171e-10  -9.78546e-11  -8.33283e-11  -7.26309e-11
 0.000849526  -2.65506e-6  -1.48303e-7  -9.15682e-8  -6.47598e-8     -1.19867e-10  -9.80994e-11  -8.21291e-11  -6.83873e-11
 0.000849248  -2.62298e-6  -1.48649e-7  -9.21537e-8  -6.55392e-8  …  -1.44725e-10  -1.07284e-10  -8.51327e-11  -7.16104e-11
 0.00084901   -2.59559e-6  -1.48948e-7  -9.26608e-8  -6.62178e-8     -2.48154e-10  -1.59995e-10  -1.08585e-10  -7.82972e-11
 0.000848811  -2.57284e-6  -1.49199e-7  -9.30869e-8  -6.67909e-8     -5.67241e-10  -3.66004e-10  -2.322e-10    -1.48171e-10
 0.000848652  -2.5547e-6   -1.494e-7    -9.34302e-8  -6.72545e-8     -1.30662e-9   -9.4078e-10   -6.57973e-10  -4.46048e-10
 0.000848533  -2.54112e-6  -1.49552e-7  -9.36891e-8  -6.76052e-8     -2.57098e-9   -2.07407e-9   -1.64632e-9   -1.28368e-9
 0.000848453  -2.53208e-6  -1.49653e-7  -9.38624e-8  -6.78404e-8  …  -4.08749e-9   -3.57702e-9   -3.11482e-9   -2.69659e-9
 0.000848413  -2.52757e-6  -1.49704e-7  -9.39492e-8  -6.79584e-8     -5.16548e-9   -4.71124e-9   -4.30053e-9   -3.92791e-9
```

```jldoctest ex
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = ContinuousWavelets.icwt(res, c, DualFrames())

# output

┌ Warning: the canonical dual frame is off by 3.8053445963886525e6, consider using one of the delta dual frames
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:37
2047-element Vector{Float64}:
  0.021353468263494175
  0.021355106707309952
  0.021358383142385133
  0.0213632966799539
  0.021369846027086553
  0.02137802953320333
  0.021387845249242986
  0.021399290997235258
  0.021412364447637153
  ⋮
 -0.006147601225469837
 -0.006147062558231826
 -0.006146592492113732
 -0.006146190563893243
 -0.006145856346051334
 -0.006145589459921974
 -0.006145389587395019
 -0.006145256480782908
 -0.006145189970537515
```

```jldoctest ex; output=false
p1 = plot(f, legend=false, title="Smoothing and dropping bumps", linewidth=2)
plot!(dropped, linewidth=3)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)

# output

Plot{Plots.PyPlotBackend() n=3}
```

![Bumps](/docs/bumps.svg)

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```jldoctest ex
julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0)


julia> res = circshift(ContinuousWavelets.cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.2677814440444114
└ @ ContinuousWavelets ~/allHail/projects/ContinuousWavelets/src/sanityChecks.jl:10
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
  1.89367e-5   0.000266033  0.000196408  2.69195e-5  -3.89652e-5   …   0.000121129   9.23682e-5    4.6727e-5     2.99983e-6
  8.33321e-5   0.000266913  0.000201712  3.29627e-5  -2.66883e-5      -1.90633e-5    3.35839e-7    1.56557e-5   -4.46419e-5
 -0.000103328  0.000267843  0.000207148  3.99906e-5  -1.42476e-5      -0.000169306  -0.000138073  -8.95268e-5   -0.000218151
 -0.000354648  0.000268836  0.000212775  4.76928e-5  -2.62643e-6      -0.000225955  -0.000244553  -0.000251492  -0.000256004
 -0.000189096  0.000269891  0.000218496  5.59397e-5   8.68511e-6      -8.89207e-5   -0.000139912  -0.000206374   5.48004e-5
  0.000341582  0.000270994  0.000224054  6.40344e-5   1.92346e-5   …   0.000160383   0.000155859   0.000132569   0.000428554
  0.000574352  0.000272142  0.000229273  7.05717e-5   2.64068e-5       0.000274473   0.000342921   0.000418284   0.000352239
  0.000103119  0.000273352  0.000234229  7.43258e-5   2.73536e-5       9.54274e-5    0.000158369   0.000244307  -0.000205831
 -0.000603783  0.000274651  0.000239099  7.51774e-5   2.16395e-5      -0.000232655  -0.000259222  -0.000284486  -0.000639977
  ⋮                                                                ⋱                               ⋮
 -4.89024e-6   0.0020357    0.00136985   0.00106421   0.000655789      2.36209e-6    6.41261e-7    1.37549e-6    1.01282e-6
 -3.99268e-6   0.00202884   0.0013532    0.00105939   0.000652748      2.63171e-6    9.60498e-7    1.47322e-6    1.13924e-6
 -3.08658e-6   0.00202194   0.00133647   0.00105447   0.000649604  …   3.00328e-6    1.37656e-6    1.66221e-6    1.36753e-6
 -2.1235e-6    0.00201502   0.00131965   0.00104945   0.000646355      3.48171e-6    1.90218e-6    1.96868e-6    1.72503e-6
 -1.08074e-6   0.00200808   0.00130275   0.00104434   0.000643004      4.06464e-6    2.53947e-6    2.40611e-6    2.2242e-6
  4.79044e-8   0.00200111   0.00128576   0.00103914   0.000639549      4.7377e-6     3.28078e-6    2.97486e-6    2.86052e-6
  1.248e-6     0.00199411   0.0012687    0.00103384   0.000635993      5.35474e-6    4.03275e-6    3.66092e-6    3.49912e-6
  2.24677e-6   0.00198709   0.00125155   0.00102844   0.000632336  …   5.66737e-6    4.55137e-6    4.24042e-6    3.80685e-6
  2.63848e-6   0.00198004   0.00123433   0.00102295   0.000628578      5.49586e-6    4.58582e-6    4.3791e-6     3.47575e-6

[:, :, 2] =
  7.81007e-18  0.0226754  0.00955729  0.00470372    0.00341201   …   3.57962e-18   9.32869e-18   3.03725e-18   3.68809e-18
 -3.47114e-18  0.022684   0.00950443  0.00469405    0.00345558      -8.67785e-19   1.73557e-18  -1.73557e-18  -3.47114e-18
 -8.8948e-18   0.0226929  0.00944925  0.00468509    0.00349665      -4.33893e-18  -8.67785e-19  -9.11175e-18  -5.64061e-18
 -3.90503e-18  0.0227021  0.00939273  0.0046771     0.00353795      -2.60336e-18  -4.33893e-18  -7.81007e-18   1.73557e-18
  5.20671e-18  0.0227115  0.00933463  0.00466972    0.00358021      -2.16946e-18   1.30168e-18  -4.77282e-18   6.72534e-18
  1.08473e-17  0.0227211  0.00927531  0.00466315    0.00362167   …   1.73557e-18   2.60336e-18   5.20671e-18   4.33893e-18
  1.73557e-18  0.0227308  0.00921533  0.00465729    0.00366155      -1.84404e-18   3.68809e-18  -2.16946e-19  -1.54574e-18
 -4.33893e-18  0.0227406  0.0091542   0.00464097    0.00369962      -2.38641e-18   2.16946e-18  -4.12198e-18  -9.11175e-18
 -6.0745e-18   0.0227506  0.00909149  0.00461836    0.00373611      -6.29144e-18  -2.16946e-18  -6.94228e-18  -3.47114e-18
  ⋮                                                              ⋱                               ⋮
 -1.16366e-17  0.0336961  0.0110364   0.000352368  -1.47154e-17     -3.58179e-18  -1.07746e-17  -1.02248e-17  -6.74285e-18
 -2.55726e-18  0.033762   0.0109992   0.000375385  -1.35714e-18     -4.72663e-19  -6.52154e-18  -3.90631e-18  -4.61867e-18
 -7.6293e-18   0.0338274  0.010964    0.000399294  -1.25056e-17  …  -2.32007e-18  -6.17772e-18  -2.21136e-18  -5.25095e-18
 -4.66963e-18  0.0338924  0.0109308   0.000424023  -8.39451e-18      3.2299e-18   -1.83791e-18   6.90007e-19  -4.72535e-18
 -3.30708e-18  0.033957   0.0108995   0.000449495  -1.65074e-17      2.13675e-18  -3.22622e-18   4.58688e-19  -7.19358e-18
 -7.6997e-18   0.0340217  0.0108677   0.000475705  -6.89573e-18      1.24614e-18  -4.3671e-18   -3.49203e-18  -6.23824e-18
 -1.03499e-17  0.0340864  0.0108359   0.000502628  -1.43084e-17      4.5178e-18   -1.82058e-18   9.12161e-19  -5.50979e-18
 -9.29595e-18  0.0341512  0.0108039   0.00053027   -1.27059e-17  …   8.83625e-19  -6.65106e-18  -3.40533e-18  -3.84208e-19
  1.27592e-18  0.0342157  0.0107729   0.000558637  -1.03715e-17      4.50903e-18  -1.1533e-18    4.21959e-18  -1.7043e-18

[:, :, 3] =
 -4.2736e-7   0.0059687   0.00256803  0.000892506  0.000666694  …  3.97814e-8   2.49033e-8   -3.38814e-8  4.47839e-8  1.86209e-8
 -4.39691e-7  0.00596762  0.00256378  0.000882834  0.000692389     2.76986e-8   1.21816e-8   -4.59854e-8  3.30771e-8  7.78201e-9
 -4.48084e-7  0.0059668   0.00256     0.000874593  0.000716669     1.61823e-8   7.4973e-10   -5.67431e-8  2.24685e-8  2.3172e-9
 -4.51634e-7  0.00596629  0.00255675  0.000868429  0.000739482     7.56909e-9  -6.89343e-9   -6.29884e-8  1.71253e-8  2.13432e-9
 -4.53736e-7  0.00596616  0.00255388  0.00086509   0.000761802     3.629e-9    -9.79742e-9   -6.45768e-8  1.66027e-8  3.47374e-9
 -4.5679e-7   0.00596647  0.00255141  0.000865675  0.00078427   …  3.56186e-9  -9.58071e-9   -6.39968e-8  1.76824e-8  4.08017e-9
 -4.6046e-7   0.0059673   0.00254963  0.000872086  0.000806868     4.59725e-9  -8.70368e-9   -6.31675e-8  1.84952e-8  4.09172e-9
 -4.64139e-7  0.00596883  0.0025489   0.000881343  0.000829423     5.69239e-9  -8.10867e-9   -6.30668e-8  1.88061e-8  4.06813e-9
 -4.67881e-7  0.00597128  0.00254952  0.000890693  0.000851916     6.29758e-9  -8.02414e-9   -6.36202e-8  1.89465e-8  4.11435e-9
  ⋮                                                             ⋱                                         ⋮
 -1.01175e-7  0.00335065  0.00169783  0.000148564  4.30591e-7      1.40432e-9  -8.62501e-10  -1.39044e-8  3.85114e-9  6.29523e-10
 -1.00662e-7  0.00335806  0.00169074  0.000151598  4.34888e-7      1.79613e-9  -5.78987e-10  -1.36661e-8  3.88963e-9  6.9526e-10
 -1.00081e-7  0.00336532  0.00168432  0.000154811  4.39299e-7   …  2.35374e-9  -1.32135e-10  -1.32681e-8  4.06365e-9  9.21894e-10
 -9.93313e-8  0.00337236  0.00167881  0.000158131  4.43829e-7      3.09208e-9   5.0037e-10   -1.26779e-8  4.42412e-9  1.36627e-9
 -9.83556e-8  0.00337911  0.00167463  0.000161542  4.48484e-7      4.03282e-9   1.34162e-9   -1.18702e-8  5.00834e-9  2.07275e-9
 -9.71063e-8  0.00338574  0.00167092  0.000165053  4.53267e-7      5.19836e-9   2.41717e-9   -1.08147e-8  5.84663e-9  3.0828e-9
 -9.55338e-8  0.00339243  0.00166694  0.000168694  4.58184e-7      6.30957e-9   3.49447e-9   -9.63413e-9  6.97634e-9  4.19423e-9
 -9.41123e-8  0.00339924  0.00166237  0.000172551  4.6324e-7    …  6.92345e-9   4.12867e-9   -8.75924e-9  8.01012e-9  4.78652e-9
 -9.36079e-8  0.0034061   0.00165751  0.000176759  4.68442e-7      6.79732e-9   3.9998e-9    -8.6509e-9   8.3188e-9   4.24252e-9

[:, :, 4] =
  0.000307454  -0.0150898   -0.00391724   0.0165289   0.0145612  …   0.000601507   0.000528795   0.000541871   0.000301757
  6.05948e-5   -0.0152536   -0.00405883   0.0162813   0.0143339      0.000347251   0.000286288   0.000307363   8.45503e-5
 -0.000106628  -0.0154172   -0.00420221   0.0160291   0.0141005      0.000118658   7.07121e-5    9.46531e-5   -2.5274e-5
 -0.000176071  -0.0155806   -0.00434733   0.0157728   0.0138617     -3.4245e-5    -5.45697e-5   -1.28804e-5   -2.95235e-5
 -0.000215585  -0.0157438   -0.00449399   0.0155126   0.0136178     -9.24075e-5   -8.66309e-5   -2.40997e-5   -3.35174e-6
 -0.000273117  -0.0159068   -0.00464219   0.0152487   0.0133693  …  -8.81642e-5   -7.52024e-5   -3.33285e-6    8.1492e-6
 -0.000341921  -0.0160696   -0.00479189   0.0149812   0.0131165     -7.06997e-5   -5.86549e-5    1.20476e-5    7.73542e-6
 -0.000409873  -0.0162322   -0.00494306   0.0147104   0.0128598     -5.88572e-5   -5.65403e-5    1.73367e-5    6.59563e-6
 -0.000477985  -0.0163945   -0.00509565   0.0144364   0.0125994     -5.72081e-5   -6.73706e-5    1.91596e-5    6.82771e-6
  ⋮                                                              ⋱                               ⋮
  0.000680261  -0.00817263  -0.00224811  -0.0184176  -0.0162279     -6.55108e-5    8.64965e-5   -2.28035e-5    7.23463e-7
  0.000611918  -0.00808368  -0.00215646  -0.018218   -0.0160653     -0.000106729   5.28847e-5   -2.86981e-5   -8.84735e-6
  0.00053323   -0.00799495  -0.0020629   -0.0180132  -0.0158957  …  -0.000171303  -3.66063e-6   -5.39545e-5   -4.14279e-5
  0.000429829  -0.00790643  -0.00196748  -0.0178031  -0.0157191     -0.000262417  -8.7768e-5    -0.000105853  -0.000105146
  0.000293499  -0.00781811  -0.00187021  -0.0175877  -0.0155353     -0.000383364  -0.000203043  -0.000189723  -0.000206338
  0.000117459  -0.00773001  -0.00177109  -0.017367   -0.0153443     -0.000537781  -0.000353818  -0.00030989   -0.000350925
 -0.000105351  -0.0076421   -0.0016701   -0.0171409  -0.0151462     -0.000692387  -0.000522514  -0.000471683  -0.000509964
 -0.000307094  -0.00755439  -0.00156729  -0.0169094  -0.0149407  …  -0.000783526  -0.000647484  -0.000619672  -0.000594673
 -0.000378125  -0.00746687  -0.00146262  -0.0166726  -0.0147279     -0.000765417  -0.000662793  -0.000663838  -0.00051676
```

And the plot of these:

```jldoctest ex; output = false
p1 = plot(plot(exs[:, 1], legend=false, title="Doppler", yticks=[], xticks=[], linewidth=2), plot(exs[:, 2], legend=false, title="Blocks", yticks=[], xticks=[], linewidth=2), plot(exs[:, 3], legend=false, title="Bumps", yticks=[], xticks=[], linewidth=2), plot(exs[:, 4], legend=false, title="HeaviSine", yticks=[], xticks=[], linewidth=2), layout=(1, 4))
p2 = plot(heatmap(identity.(res[:, :, 1])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 2])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 3])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 4])', xticks=false, yticks=[], c=:viridis, colorbar=false), layout=(1, 4))
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)
```

![parallel transforms](/docs/multiEx.svg)

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
