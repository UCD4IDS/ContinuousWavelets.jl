
function testFourierDomainProperties(daughters, isAve)
    dMags = abs.(daughters)
    lowAprxAnalyt = daughters[1, 2] ./ maximum(dMags[:, 2])
    if abs(lowAprxAnalyt) >= 0.01
        lowAprxAnalyt = round(lowAprxAnalyt, sigdigits=5)
        @warn "the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully" lowAprxAnalyt
    end
    highAprxAnalyt = abs.(daughters[end, end]) ./ maximum(dMags[:, 2])
    if highAprxAnalyt >= 0.01
        highAprxAnalyt = round(highAprxAnalyt, sigdigits=5)
        @warn "the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully" highAprxAnalyt
    end

    netWeight = sum(dMags, dims = 2)
    peakArgs = argmax(dMags[:, (isAve+1):end], dims = 1)
    peakFreqs = dMags[:, isAve+1:end][peakArgs]
    maxPeak = maximum(peakFreqs)
    minPeak = minimum(peakFreqs)
    ratioOfCoverage = maxPeak / minPeak
    if ratioOfCoverage > 5
        ratioOfCoverage = round(ratioOfCoverage, sigdigits=5)
        @warn "there are frequencies which are significantly more covered than others, with a ratio of" ratioOfCoverage
    end
    peakLocs = [x[1] for x in peakArgs] # the location of the peaks, irrespective of which wavelet has that peak location
    first = minimum(peakLocs) # the first peak
    last = maximum(peakLocs)
    minimalRegionComparedToLastPeak = maximum(peakFreqs) /
                                      minimum(abs.(netWeight[first:last]))
    if minimalRegionComparedToLastPeak > 2
        minimalRegionComparedToLastPeak = round(minimalRegionComparedToLastPeak, sigdigits=5)
        @warn "there are wavelets whose peaks are far enough apart that the trough between them is less than half the height of the highest frequency wavelet" minimalRegionComparedToLastPeak
    end
end


function testDualCoverage(canonicalDualFrame)
    dualCover = sum(canonicalDualFrame, dims = 2)
    err = norm(dualCover .- 1)
    if err > 50
        @warn "the canonical dual frame is off by $(err), consider using one of the delta dual frames"
    end
end

function testDualCoverage(β, Ŵ)
    dualCover = sum(conj.(β) .* Ŵ, dims = 2)
    err = norm(dualCover .- 1)
    if err > 50
        @warn "the dual frame is off by $(err), consider using a different dual frame"
    end
end
