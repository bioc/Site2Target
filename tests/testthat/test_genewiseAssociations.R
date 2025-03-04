suppressPackageStartupMessages({
    library(Site2Target)
})

test_that("Remove reserve character test",{
    expect_identical(removeReserveCharacter("A&%B^f6"),"ABf6")
})
