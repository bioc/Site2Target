suppressPackageStartupMessages({
    library(Site2Target)
})

test_that("converting-Centering string to granges and vice versa",{
    peakStr <- c("chr1:100-130", "chr2:200-220")
    peakStrCentered <- c("chr1:115-115", "chr2:210-210")
    peakGR <- string2Granges(peakStr)
    expect_identical(granges2String(peakGR),peakStr)
    peakGRCentered <- getCenterOfPeaks(peakGR)
    expect_identical(granges2String(peakGRCentered),peakStrCentered)
})

