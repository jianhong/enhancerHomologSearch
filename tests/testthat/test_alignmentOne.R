test_that("alignmentOne works not correct", {
  suppressMessages(peaks <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  peaks <- peaks[seqnames(peaks)=="chrX" & width(peaks)<10000]
  peaks$id <- paste(names(peaks), 1, sep="_")
  subj <- Enhancers(genome=Hsapiens, peaks=peaks)
  q <- getSeq(Hsapiens, peaks[which.min(width(peaks))]) # MIR320D2
  bpparam <- MulticoreParam(workers = 1, tasks=200, progressbar=TRUE)
  ao <- alignmentOne(q, subj, bpparam=bpparam)
  res <- peaks(ao)
  res <- res[which.max(res$score)]
  expect_equal(res$gene_id, names(q))
})
