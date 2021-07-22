test_that("guessAssembly works not correct", {
    assembly <- enhancerHomologSearch:::guessAssembly(Hsapiens)
    expect_is(assembly, "character")
    expect_equal(assembly, c("GRCh38", "hg38"))

    assembly <- enhancerHomologSearch:::guessAssembly(Mmusculus)
    expect_equal(assembly, c("GRCm38", "mm10"))
})
