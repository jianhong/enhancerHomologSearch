#' Transcription Factor Binding Pattern Similarity (TFBPS) search
#' @description Search the TFBPs for query in subject.
#' @param query An object of DNAStringSet to represent enhancer
#' @param subject Output of getENCODEdata. An object of \link{Enhancers}
#' @param PWMs The Position Weight Matrix list represented as a numeric matrix.
#' Object of \link[TFBSTools:XMatrixList]{PWMatrixList} or
#' \link[TFBSTools:XMatrixList]{PFMatrixList}.
#' @param queryGenome An object of \link[BSgenome:BSgenome-class]{BSgenome} for
#' query data.
#' @param background background nucleotide frequencies. Default is "genome".
#' Refer \link[motifmatchr]{matchMotifs} for details.
#' @param \dots Parameters will be passed to \link[motifmatchr]{matchMotifs}
#' except 'out' and 'genome'.
#' @param maximalShuffleEnhancers The maximal number of Shuffled enhancers.
#' If the number of the input enhancer candidates is greater than
#' maximalShuffleEnhancers, no shuffled enhancer sequences will be included.
#' The shuffled enhancers will be created by \link{shuffle}.
#' @return An object of \link{Enhancers}.
#' @importFrom BiocGenerics score
#' @importFrom IRanges subject
#' @importFrom stats p.adjust
#' @export
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' peaks <- GRanges("chr1", IRanges(seq(5000, 50000, by=1000), width=1000))
#' peaks$id <- paste(seq_along(peaks), 1, sep="_")
#' subj <- Enhancers(genome=Hsapiens, peaks=peaks)
#' q <- getSeq(Hsapiens, GRanges("chr1", IRanges(90000, width=1000)))
#' data(motifs)
#' ao <- searchTFBPS(q, subj, motifs[["dist60"]], queryGenome=Hsapiens,
#'                   maximalShuffleEnhancers = 50)
searchTFBPS <- function(query, subject, PWMs, queryGenome,
                        background="genome", ...,
                        maximalShuffleEnhancers=1000){
  checkQuerySubject(query, subject, subjectIsList=FALSE)
  peaks <- subject@peaks
  checkPWMs(PWMs)
  background <- match.arg(background, choices = c("subject", "genome", "even"))
  ## search TFBP in query
  target <- getSeq(subject, peaks)
  if(length(target)<maximalShuffleEnhancers){
    shuffled <- shuffle(target[grepl("^fwd", names(target))],
                        n=ceiling(2*maximalShuffleEnhancers/length(target)))
    shuffled <- sample(shuffled,
                       size = maximalShuffleEnhancers - length(target))
    target <- c(target, shuffled)
  }

  TFBP_query <- searchTFBP(query, PWMs,
                           genome=queryGenome,
                           bg=background, ...)[1, ]
  TFBP_subject <- searchTFBP(target, PWMs,
                             genome=genome(subject),
                             bg=background, ...)
  qid <- rownames(TFBP_subject)
  qid[!grepl("shuffle", qid)] <-
    paste0(qid[!grepl("shuffle", qid)], "_target_0")
  qid <- do.call(rbind, strsplit(qid, "_"))
  stopifnot("subject must be output of getENCODEdata."=ncol(qid)==5)
  cnqid <- c("dire", "peakID", "tileID", "type", "idx")
  colnames(qid) <- cnqid
  scores <- apply(TFBP_subject, 1, function(.ele){
    TFBPscore(TFBP_query, .ele)
  })
  qid <- as.data.frame(qid)
  qid <- cbind(qid, score=scores)
  qid$order <- seq.int(nrow(qid))
  qid <- qid[order(qid$score, decreasing = TRUE), , drop=FALSE]
  qid <- qid[!duplicated(paste(qid$peakID, qid$type, qid$idx)), , drop=FALSE]
  qid <- qid[order(qid$order), , drop=FALSE]
  ## get Z-score
  m <- mean(qid$score)
  sd <- sd(qid$score)
  qid$Z <- (qid$score - m)/sd

  S <- nrow(qid)
  qid$pval <- vapply(qid$score, FUN=function(.ele){
    sum(qid$score>=.ele)/S
  }, FUN.VALUE = numeric(1))
  qid$adjp <- p.adjust(qid$pval, method = "BH")
  qid <- qid[qid$type=="target", , drop=FALSE]

  peaks <- peaks[match(paste(qid$peakID, qid$tileID, sep="_"), peaks$id)]
  strand(peaks) <- ifelse(qid$dire=="fwd",
                           ifelse(strand(peaks)=="-", "-", "+"),
                           ifelse(strand(peaks)=="-", "+", "-"))
  for(cn in c("score", "Z", "pval", "adjp")){
    mcols(peaks)[[cn]] <- qid[[cn]]
  }
  TFBP <-TFBP_subject[apply(qid[, cnqid[seq.int(3)]], 1,
                            paste, collapse = "_"),
                       , drop=FALSE]
  rownames(TFBP) <- sub("(fwd|rev)_", "", rownames(TFBP))

  Enhancers(genome = subject@genome, peaks = peaks, TFBP=TFBP, TFBP0=TFBP_query)
}
## input DNAStringSet and PWMs
## output PWM match pattern: barcode of match.
## the barcode is a named logical vector
#' @importFrom motifmatchr matchMotifs motifMatches
searchTFBP <- function(x, pwms, ...){
  stopifnot(is(x, "DNAStringSet"))
  checkPWMs(pwms)
  ## barcode for data
  args <- list(subject=x, pwms=pwms, ...)
  args$out <- "matches"
  TFBP <- do.call(matchMotifs, args = args)
  mm <- motifMatches(TFBP)
  rownames(mm) <- names(x)
  args_rev <- args
  args_rev$subject <- reverseComplement(x)
  TFBP_rev <- do.call(matchMotifs, args = args_rev)
  mm_rev <- motifMatches(TFBP_rev)
  rownames(mm_rev) <- names(x)
  mm <- mm | mm_rev
  return(mm)
}
## compare TFBP between query and subjects
TFBPscore <- function(query, subject, weight=NULL){
  checkTFBPs(query)
  checkTFBPs(subject)
  share <- intersect(names(query), names(subject))
  equal <- query[share] & subject[share]
  alln <- union(names(query), names(subject))
  if(all(alln %in% names(weight))){
    equal <- equal*weight[share]
    alln <- sum(weight[alln])
  }else{
    alln <- length(alln)
  }
  JaccardIndex <- sum(equal)/alln
  return(JaccardIndex)
}
