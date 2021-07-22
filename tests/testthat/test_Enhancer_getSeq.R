test_that("Enhancer::getSeq works not correct", {
  suppressMessages(peaks <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  enh <- Enhancers(genome=Hsapiens, peaks=peaks[1:5])
  seq1 <- getSeq(Hsapiens, peaks[1:5])
  seq2 <- getSeq(enh)
  expect_equal(unname(c(seq1, reverseComplement(seq1))),
               unname(seq2))

  seq <- getSeq(enh, peaks[1])
  expect_equal(length(seq), 2)
  expect_true(nchar(as.character(seq))[1]==width(peaks)[1])
})
