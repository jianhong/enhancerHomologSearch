library(MotifDb)
library(motifStack)
library(TFBSTools)
library(ade4)
motifs <- query(MotifDb, "hsapiens|mmusculus")
MotifList2PWMatrixList <- function(motifs, groupDistance=10){
  motifs <- motifs[!is.na(mcols(motifs)$geneSymbol)]
  meta <- mcols(motifs)
  mat <- as.list(motifs)
  keep <- meta$dataSource %in% c("jaspar2018", "jolma2013", "cisbp_1.02")
  meta <- meta[keep, ]
  mat <- mat[keep]
  oid <- order(c("jaspar2018"=2, "jolma2013"=1, "cisbp_1.02"=3)[meta$dataSource])
  meta <- meta[oid, ]
  mat <- mat[oid]
  keep <- !duplicated(tolower(meta$geneSymbol))
  meta <- meta[keep, ]
  mat <- mat[keep]
  align <- matalign(mat)
  hc <- motifHclust(align, method="average")
  phylog <- hclust2phylog(hc)
  leaves <- names(phylog$leaves)
  names(mat) <- gsub("[^a-zA-Z0-9]","_",
                     gsub("(_[0-9]+)+$", "", names(mat)))
  rownames(meta) <- gsub("[^a-zA-Z0-9]","_",
                         gsub("(_[0-9]+)+$", "", rownames(meta)))
  pfms <- mat[leaves]
  pfms <- lapply(names(pfms), function(.ele, pfms){new("pfm",mat=pfms[[.ele]],
                                                       name=.ele)},pfms)
  out <- lapply(groupDistance, function(gD){
    motifSig <- motifSignature(pfms, phylog,
                               groupDistance=gD,
                               min.freq = 1)
    # pdf(paste0("motifPiles.", gD, ".pdf"), height = 40, width=12)
    # motifPiles(phylog=phylog, pfms=DNAmotifAlignment(pfms),
    #            pfms2=signatures(motifSig),
    #            motifScale="logarithmic",
    #            col.pfms2=sigColor(motifSig),
    #            col.pfms2.width=.01,
    #            plotIndex=TRUE,
    #            groupDistance=gD)
    # dev.off()
    sig <- signatures(motifSig)
    newMotifs <- lapply(sig, function(.ele){
      n <- .ele@name
      n <- strsplit(n, ";")[[1]][1]
      m <- mat[[n]]
      m[4, ] <- 1 - colSums(m[-4, ])
      tag <- meta[n, , drop=TRUE]
      PWMatrix(ID=n, name=tag$geneSymbol, matrixClass=tag$tfFamily, tags=tag,
               profileMatrix=m)
    })
    names(newMotifs) <- sapply(newMotifs, function(.ele) .ele@name)
    do.call(PWMatrixList, newMotifs)
  })
  return(out)
}
motifs <- MotifList2PWMatrixList(motifs, seq(10, 100, by=10))
names(motifs) <- paste0("dist", seq(10, 100, by=10))
saveRDS(motifs, "PWMatrixList.rds")
