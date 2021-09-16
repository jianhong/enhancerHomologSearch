## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
suppressPackageStartupMessages({
  library(enhancerHomologSearch)
  library(BSgenome.Drerio.UCSC.danRer10)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

## -----------------------------------------------------------------------------
# load genome sequences
library(BSgenome.Drerio.UCSC.danRer10)
# define the enhancer genomic coordinates
LEN <- GRanges("chr4", IRanges(19050041, 19051709))
# extract the sequences as Biostrings::DNAStringSet object
(seqEN <- getSeq(BSgenome.Drerio.UCSC.danRer10, LEN))

## -----------------------------------------------------------------------------
# load library
library(enhancerHomologSearch)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
# download enhancer candidates for human heart tissue
hs <- getENCODEdata(genome=Hsapiens,
                    partialMatch=c(biosample_summary = "heart"))
# download enhancer candidates for mouse heart tissue
mm <- getENCODEdata(genome=Mmusculus,
                    partialMatch=c(biosample_summary = "heart"))

## -----------------------------------------------------------------------------
# subset the data for test run
# In this test run, we will only use upstream 1M and downstream 1M of homolog
# gene
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
eid <- mget("LEP", org.Hs.egALIAS2EG)[[1]]
g_hs <- select(TxDb.Hsapiens.UCSC.hg38.knownGene,
               keys=eid,
               columns=c("GENEID", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND"),
               keytype="GENEID")
g_hs <- range(with(g_hs, GRanges(TXCHROM, IRanges(TXSTART, TXEND))))
expandGR <- function(x, ext){
  stopifnot(length(x)==1)
  start(x) <- max(1, start(x)-ext)
  end(x) <- end(x)+ext
  GenomicRanges::trim(x)
}
hs <- subsetByOverlaps(hs, expandGR(g_hs, ext=1000000))
# Here we use the subset of 1M upstream and downstream of homolog gene.
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
eid <- mget("Lep", org.Mm.egALIAS2EG)[[1]]
g_mm <- select(TxDb.Mmusculus.UCSC.mm10.knownGene,
               keys=eid,
               columns=c("GENEID", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND"),
               keytype="GENEID")
g_mm <- range(with(g_mm,
                   GRanges(TXCHROM,
                           IRanges(TXSTART, TXEND),
                           strand=TXSTRAND)))
g_mm <- g_mm[seqnames(g_mm) %in% "chr6" & strand(g_mm) %in% "+"]
mm <- subsetByOverlaps(mm, expandGR(g_mm, ext=1000000))

# search the binding pattern
data(motifs)
PWMs <- motifs[["dist60"]]
aln_hs <- searchTFBPS(seqEN, hs, PWMs = PWMs, queryGenome = Drerio)
aln_mm <- searchTFBPS(seqEN, mm, PWMs = PWMs, queryGenome = Drerio)

## -----------------------------------------------------------------------------
# Step4
ext <- 100000
aln_hs <- subsetByOverlaps(aln_hs, ranges = expandGR(g_hs, ext=ext))
## filter by distance
distance(aln_hs) <- distance(peaks(aln_hs), g_hs, ignore.strand=TRUE)
subset(aln_hs, pval<0.1 & distance >5000)
aln_hs <- subset(aln_hs, pval<0.1 & distance >5000)
saveRDS(aln_hs, "inst/extdata/aln_hs.rds")


aln_mm <- subsetByOverlaps(aln_mm, ranges = expandGR(g_mm, ext=ext))
## filter by distance
distance(aln_mm) <- distance(peaks(aln_mm), g_mm, ignore.strand=TRUE)
subset(aln_mm, pval<0.1 & distance >5000)
aln_mm <- subset(aln_mm, pval<0.1 & distance >5000)
saveRDS(aln_mm, "inst/extdata/aln_mm.rds")

## -----------------------------------------------------------------------------
al <- alignment(seqEN, list(human=aln_hs, mouse=aln_mm),
                method="ClustalW", order="input")
saveRDS(al, "inst/extdata/al.rds")
