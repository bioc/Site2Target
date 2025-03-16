suppressPackageStartupMessages({
    library(Site2Target)
})

geneFile=system.file("extdata", "gene_expression.tsv",
                     package="Site2Target")
geneCoords <- Table2Granges(geneFile)
geneTable <- read.table(geneFile, header=TRUE)

tfFile =system.file("extdata", "MEIS_binding.tsv",
                    package="Site2Target")
TFCoords <- Table2Granges(tfFile)
tfTable <- read.table(tfFile, header=TRUE)
tfIntensities <- tfTable$intensities

HiCFile =system.file("extdata", "HiC_intensities.tsv",
                     package="Site2Target")
HiCstr1 <- Table2Granges(HiCFile, chrColName="Strand1_chr",
                         startColName="Strand1_start", endColName="Strand1_end")
HiCstr2 <- Table2Granges(HiCFile, chrColName="Strand2_chr",
                         startColName="Strand2_start", endColName="Strand2_end")


test_that("Check getTargetGenesPvals",{
    pvals <- getTargetGenesPvals( geneCoordinates=geneCoords,
                                  sites=TFCoords, distance = 50000)

    topTargetIndex <- order(pvals)[1]
    expect_identical(geneTable$name[topTargetIndex], "RUNX1")
    expect_equal(pvals[topTargetIndex], 2.22594825990949e-06, tolerance = 1e-6)

})

test_that("Check getTargetGenesPvalsWithDNAInteractions",{
    pvals <-
        getTargetGenesPvalsWithDNAInteractions(
            geneCoordinates=geneCoords,
            sites=TFCoords,
            strand1=HiCstr1,
            strand2=HiCstr2)

    topTargetIndex <- order(pvals)[1]
    expect_identical(geneTable$name[topTargetIndex], "RUNX1")
    expect_equal(pvals[topTargetIndex], 2.22594825990949e-06, tolerance = 1e-6)

})

test_that("Check getTargetGenesPvalsWithIntensities",{
    pvals <-
        getTargetGenesPvalsWithIntensities(
            geneCoordinates=geneCoords,
            sites=TFCoords,
            intensities=tfIntensities
        )

    topTargetIndex <- order(pvals)[1]
    expect_identical(geneTable$name[topTargetIndex], "RUNX1")
    expect_equal(pvals[topTargetIndex], 0.00169937052797075, tolerance = 1e-4)

})


test_that("Check getTargetGenesPvalsWithIntensitiesAndDNAInteractions",{
    pvals <-
        getTargetGenesPvalsWithIntensitiesAndDNAInteractions(
            geneCoordinates=geneCoords,
            sites=TFCoords,
            intensities=tfIntensities,
            strand1=HiCstr1,
            strand2=HiCstr2
        )

    topTargetIndex <- order(pvals)[1]
    expect_identical(geneTable$name[topTargetIndex], "RUNX1")
    expect_equal(pvals[topTargetIndex], 0.00169937052797075, tolerance = 1e-4)

})
