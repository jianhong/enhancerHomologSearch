test_that("searchTFBPS works not correct", {
  suppressMessages(peaks <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  peaks <- peaks[seqnames(peaks)=="chrX" & width(peaks)<3000 &
                   width(peaks)>500]
  peaks$id <- paste(names(peaks), 1, sep="_")
  subj <- Enhancers(genome=Hsapiens, peaks=peaks)
  q <- getSeq(Hsapiens, peaks[which.min(width(peaks))]) # MIR320D2
  data(motifs)
  ao <- searchTFBPS(q, subj, PWMs=motifs[["dist60"]],
                    queryGenome = Hsapiens)
  res <- peaks(ao)
  expect_equal(res[names(q)]$score, max(res$score))
})
