```@meta
DocTestFilters = r"\@ ContinuousWavelets .*"
```

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

```jldoctest ex; filter= r"[0-9]\.[0-9]{5}e-[0-9][5-9]"
julia> using Random

julia> Random.seed!(1234);

julia> using ContinuousWavelets, Wavelets

julia> n = 2047;

julia> t = range(0, n / 1000, length=n); # 1kHz sampling rate

julia> f = testfunction(n, "Doppler");

julia> c = wavelet(Morlet(π), β=2)

julia> res = ContinuousWavelets.cwt(f, c)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:6
2047×31 Matrix{ComplexF64}:
 -2.84943e-6+3.86536e-19im  …   0.000109884+9.67013e-5im
 -2.84699e-6-6.11361e-20im      -8.24222e-5+0.000130545im
 -2.84212e-6+4.37411e-20im     -0.000153333-5.64666e-5im
 -2.83483e-6-3.11387e-19im       1.90839e-5-0.00016841im
 -2.82514e-6-1.31096e-19im      0.000172696-2.56466e-5im
 -2.81306e-6-3.38731e-19im  …    7.79501e-5+0.000162848im
 -2.79865e-6-9.8753e-19im      -0.000128919+0.000132755im
 -2.78192e-6+4.91292e-20im     -0.000172323-6.71036e-5im
 -2.76293e-6+5.80924e-19im      -9.39619e-6-0.000179998im
 -2.74172e-6+1.11752e-19im      0.000153988-8.14305e-5im
            ⋮               ⋱              ⋮
 0.000172941+2.7571e-19im       -2.59966e-6-7.0039e-7im
 0.000171274+1.41585e-19im      -2.58362e-6-6.24181e-7im
 0.000169814-7.90531e-21im  …   -2.56497e-6-5.41108e-7im
 0.000168561-5.81895e-20im      -2.56135e-6-4.4822e-7im
 0.000167516-3.07438e-19im      -2.57801e-6-3.83539e-7im
 0.000166679-6.64104e-19im      -2.54192e-6-3.54776e-7im
 0.000166051-1.45091e-18im       -2.4431e-6-2.37279e-7im
 0.000165633+5.67188e-19im  …   -2.46986e-6-1.96881e-8im
 0.000165423+1.25225e-18im      -2.63276e-6+4.61939e-8im


```

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c)

# output

2047×27 Matrix{Float64}:
 0.00262067   -0.00150482  -2.16134e-6  …  -3.11967e-8   -2.84659e-8
 0.00262075   -0.00150528  -2.16503e-6     -2.27191e-8   -1.96643e-8
 0.00262089   -0.00150621  -2.17246e-6     -1.22211e-8   -9.55694e-9
 0.00262111   -0.0015076   -2.18368e-6     -5.16474e-9   -3.5861e-9
 0.00262141   -0.00150945  -2.1988e-6      -2.13149e-9   -1.45404e-9
 0.00262177   -0.00151176  -2.21794e-6  …  -1.26651e-9   -9.85913e-10
 0.00262221   -0.00151453  -2.24129e-6     -1.11251e-9   -9.3166e-10
 0.00262272   -0.00151775  -2.26904e-6     -1.11305e-9   -9.48609e-10
 0.00262331   -0.00152143  -2.30144e-6     -1.13556e-9   -9.66579e-10
 0.00262397   -0.00152557  -2.33879e-6     -1.1622e-9    -9.92573e-10
 ⋮                                      ⋱   ⋮
 0.000849844  -2.6919e-6   -1.4791e-7      -8.33283e-11  -7.26309e-11
 0.000849526  -2.65506e-6  -1.48303e-7     -8.21291e-11  -6.83873e-11
 0.000849248  -2.62298e-6  -1.48649e-7  …  -8.51327e-11  -7.16104e-11
 0.00084901   -2.59559e-6  -1.48948e-7     -1.08585e-10  -7.82972e-11
 0.000848811  -2.57284e-6  -1.49199e-7     -2.322e-10    -1.48171e-10
 0.000848652  -2.5547e-6   -1.494e-7       -6.57973e-10  -4.46048e-10
 0.000848533  -2.54112e-6  -1.49552e-7     -1.64632e-9   -1.28368e-9
 0.000848453  -2.53208e-6  -1.49653e-7  …  -3.11482e-9   -2.69659e-9
 0.000848413  -2.52757e-6  -1.49704e-7     -4.30053e-9   -3.92791e-9
```

```jldoctest ex
using ContinuousWavelets, Wavelets
f = testfunction(n, "Bumps");
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
dropped = ContinuousWavelets.icwt(res, c, DualFrames())

# output

┌ Warning: the canonical dual frame is off by 3.8053445963886525e6, consider using one of the delta dual frames
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/src/sanityChecks.jl:37
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
  0.021427063201564893
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

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```jldoctest ex; filter= r"[0-9]\.[0-9]{5}e-[0-9][5-9]"
julia> using Wavelets

julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0)

julia> res = circshift(ContinuousWavelets.cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.2677814440444114
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:10
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
  1.89367e-5   0.000266033  0.000196408  …   4.6727e-5     2.99983e-6
  8.33321e-5   0.000266913  0.000201712      1.56557e-5   -4.46419e-5
 -0.000103328  0.000267843  0.000207148     -8.95268e-5   -0.000218151
 -0.000354648  0.000268836  0.000212775     -0.000251492  -0.000256004
 -0.000189096  0.000269891  0.000218496     -0.000206374   5.48004e-5
  0.000341582  0.000270994  0.000224054  …   0.000132569   0.000428554
  0.000574352  0.000272142  0.000229273      0.000418284   0.000352239
  0.000103119  0.000273352  0.000234229      0.000244307  -0.000205831
 -0.000603783  0.000274651  0.000239099     -0.000284486  -0.000639977
 -0.000707411  0.000276047  0.000243868     -0.000581431  -0.000390285
  ⋮                                      ⋱   ⋮
 -4.89024e-6   0.0020357    0.00136985       1.37549e-6    1.01282e-6
 -3.99268e-6   0.00202884   0.0013532        1.47322e-6    1.13924e-6
 -3.08658e-6   0.00202194   0.00133647   …   1.66221e-6    1.36753e-6
 -2.1235e-6    0.00201502   0.00131965       1.96868e-6    1.72503e-6
 -1.08074e-6   0.00200808   0.00130275       2.40611e-6    2.2242e-6
  4.79044e-8   0.00200111   0.00128576       2.97486e-6    2.86052e-6
  1.248e-6     0.00199411   0.0012687        3.66092e-6    3.49912e-6
  2.24677e-6   0.00198709   0.00125155   …   4.24042e-6    3.80685e-6
  2.63848e-6   0.00198004   0.00123433       4.3791e-6     3.47575e-6

[:, :, 2] =
  7.81007e-18  0.0226754  0.00955729  …   3.03725e-18   3.68809e-18
 -3.47114e-18  0.022684   0.00950443     -1.73557e-18  -3.47114e-18
 -8.8948e-18   0.0226929  0.00944925     -9.11175e-18  -5.64061e-18
 -3.90503e-18  0.0227021  0.00939273     -7.81007e-18   1.73557e-18
  5.20671e-18  0.0227115  0.00933463     -4.77282e-18   6.72534e-18
  1.08473e-17  0.0227211  0.00927531  …   5.20671e-18   4.33893e-18
  1.73557e-18  0.0227308  0.00921533     -2.16946e-19  -1.54574e-18
 -4.33893e-18  0.0227406  0.0091542      -4.12198e-18  -9.11175e-18
 -6.0745e-18   0.0227506  0.00909149     -6.94228e-18  -3.47114e-18
  1.73557e-18  0.0227606  0.0090268      -1.95252e-18   1.51862e-18
  ⋮                                   ⋱   ⋮
 -1.16366e-17  0.0336961  0.0110364      -1.02248e-17  -6.74285e-18
 -2.55726e-18  0.033762   0.0109992      -3.90631e-18  -4.61867e-18
 -7.6293e-18   0.0338274  0.010964    …  -2.21136e-18  -5.25095e-18
 -4.66963e-18  0.0338924  0.0109308       6.90007e-19  -4.72535e-18
 -3.30708e-18  0.033957   0.0108995       4.58688e-19  -7.19358e-18
 -7.6997e-18   0.0340217  0.0108677      -3.49203e-18  -6.23824e-18
 -1.03499e-17  0.0340864  0.0108359       9.12161e-19  -5.50979e-18
 -9.29595e-18  0.0341512  0.0108039   …  -3.40533e-18  -3.84208e-19
  1.27592e-18  0.0342157  0.0107729       4.21959e-18  -1.7043e-18

[:, :, 3] =
 -4.2736e-7   0.0059687   0.00256803  0.000892506  …  4.47839e-8  1.86209e-8
 -4.39691e-7  0.00596762  0.00256378  0.000882834     3.30771e-8  7.78201e-9
 -4.48084e-7  0.0059668   0.00256     0.000874593     2.24685e-8  2.3172e-9
 -4.51634e-7  0.00596629  0.00255675  0.000868429     1.71253e-8  2.13432e-9
 -4.53736e-7  0.00596616  0.00255388  0.00086509      1.66027e-8  3.47374e-9
 -4.5679e-7   0.00596647  0.00255141  0.000865675  …  1.76824e-8  4.08017e-9
 -4.6046e-7   0.0059673   0.00254963  0.000872086     1.84952e-8  4.09172e-9
 -4.64139e-7  0.00596883  0.0025489   0.000881343     1.88061e-8  4.06813e-9
 -4.67881e-7  0.00597128  0.00254952  0.000890693     1.89465e-8  4.11435e-9
 -4.71676e-7  0.00597491  0.00255138  0.000899523     1.91417e-8  4.16722e-9
  ⋮                                                ⋱  ⋮
 -1.01175e-7  0.00335065  0.00169783  0.000148564     3.85114e-9  6.29523e-10
 -1.00662e-7  0.00335806  0.00169074  0.000151598     3.88963e-9  6.9526e-10
 -1.00081e-7  0.00336532  0.00168432  0.000154811  …  4.06365e-9  9.21894e-10
 -9.93313e-8  0.00337236  0.00167881  0.000158131     4.42412e-9  1.36627e-9
 -9.83556e-8  0.00337911  0.00167463  0.000161542     5.00834e-9  2.07275e-9
 -9.71063e-8  0.00338574  0.00167092  0.000165053     5.84663e-9  3.0828e-9
 -9.55338e-8  0.00339243  0.00166694  0.000168694     6.97634e-9  4.19423e-9
 -9.41123e-8  0.00339924  0.00166237  0.000172551  …  8.01012e-9  4.78652e-9
 -9.36079e-8  0.0034061   0.00165751  0.000176759     8.3188e-9   4.24252e-9

[:, :, 4] =
  0.000307454  -0.0150898   -0.00391724  …   0.000541871   0.000301757
  6.05948e-5   -0.0152536   -0.00405883      0.000307363   8.45503e-5
 -0.000106628  -0.0154172   -0.00420221      9.46531e-5   -2.5274e-5
 -0.000176071  -0.0155806   -0.00434733     -1.28804e-5   -2.95235e-5
 -0.000215585  -0.0157438   -0.00449399     -2.40997e-5   -3.35174e-6
 -0.000273117  -0.0159068   -0.00464219  …  -3.33285e-6    8.1492e-6
 -0.000341921  -0.0160696   -0.00479189      1.20476e-5    7.73542e-6
 -0.000409873  -0.0162322   -0.00494306      1.73367e-5    6.59563e-6
 -0.000477985  -0.0163945   -0.00509565      1.91596e-5    6.82771e-6
 -0.000546078  -0.0165567   -0.00524962      2.20147e-5    7.16823e-6
  ⋮                                      ⋱   ⋮
  0.000680261  -0.00817263  -0.00224811     -2.28035e-5    7.23463e-7
  0.000611918  -0.00808368  -0.00215646     -2.86981e-5   -8.84735e-6
  0.00053323   -0.00799495  -0.0020629   …  -5.39545e-5   -4.14279e-5
  0.000429829  -0.00790643  -0.00196748     -0.000105853  -0.000105146
  0.000293499  -0.00781811  -0.00187021     -0.000189723  -0.000206338
  0.000117459  -0.00773001  -0.00177109     -0.00030989   -0.000350925
 -0.000105351  -0.0076421   -0.0016701      -0.000471683  -0.000509964
 -0.000307094  -0.00755439  -0.00156729  …  -0.000619672  -0.000594673
 -0.000378125  -0.00746687  -0.00146262     -0.000663838  -0.00051676
```

And the plot of these:

There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
