
function testFourierDomainProperties(daughters, isAve)
    lowAprxAnalyt = daughters[1,2]./maximum(abs.(daughters[:,2]))
    if abs(lowAprxAnalyt) >= .01
        @warn "the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully" lowAprxAnalyt
    end
    highAprxAnalyt = abs.(daughters[end,end])./maximum(abs.(daughters[:,2]))
    if highAprxAnalyt >= .01
        @warn "the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully" highAprxAnalyt
    end

    netWeight = sum(daughters, dims=2)
    centralFreqLast = argmax(abs.(daughters[:,end]))
    centralFreqFirst = argmax(abs.(daughters[:,1+isAve]))
    ratioOfCoverage = maximum(abs.(netWeight)) /
        minimum(abs.(netWeight[centralFreqFirst:centralFreqLast]))
    if ratioOfCoverage > 5
        @warn "there are frequencies which are significantly more covered than others, with a ratio of" ratioOfCoverage
    end
    minimalRegionComparedToLastPeak = maximum(abs.(daughters[:,end])) /
        minimum(abs.(netWeight[centralFreqFirst:centralFreqLast]))
    if minimalRegionComparedToLastPeak > 2
        @warn "there are wavelets whose peaks are far enough apart that the trough between them is less than half the height of the highest frequency wavelet" minimalRegionComparedToLastPeak
    end
end
