#' Get alignment scores
#' @description Do pairwise alignment for query enhancer to target genome
#' @param query An object of DNAStringSet to represent enhancer
#' @param subject Output of getENCODEdata. An object of \link{enhancers}
#' @param block The size of sequences to do alignment. Increase the size will
#' increase the memory cost. Default 1000.
#' @param bpparam BiocParallel parameters.
#' @param \dots not used.
#' @return An object of \link{enhancers}.
#' @importFrom Biostrings pairwiseAlignment reverseComplement unaligned
#' @importFrom BiocGenerics score
#' @importFrom IRanges subject
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
#' @export
#' @examples
#' library(BiocParallel)
#' bpparam <- MulticoreParam(workers = 2, tasks=200, progressbar=TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' peaks <- GRanges("chr1", IRanges(seq(5000, 50000, by=1000), width=1000))
#' peaks$id <- paste(seq_along(peaks), 1, sep="_")
#' subj <- enhancers(genome=Hsapiens, peaks=peaks)
#' q <- getSeq(Hsapiens, GRanges("chr1", IRanges(90000, width=1000)))
#' ao <- alignmentOne(q, subj, bpparam=bpparam)
alignmentOne <- function(query, subject, block=1000, bpparam = bpparam(), ...){
  stopifnot("query must be an object of DNAStringSet" =
              is(query, "DNAStringSet"))
  stopifnot("The length of query must be 1" = length(query)==1)
  stopifnot("subject must be an object of enhancers" =
              is(subject, "enhancers"))
  peaks <- subject@peaks
  l <- length(peaks)
  pid <- rep(seq.int(ceiling(l/block)), each = block)[seq.int(l)]
  peaks <- split(peaks, pid)
  peaks <- bplapply(X=peaks, FUN = function(.peaks, subject, query){
    target <- getSeq(subject, .peaks)
    al <- pairwiseAlignment(rep(query, length(target)),
                            subject = target, type = "global")
    query_neg <- reverseComplement(query)
    al_neg <- pairwiseAlignment(rep(query_neg, length(target)),
                                subject = target, type = "global")
    ## select the best match for each candidates region
    qscore <- ifelse(score(al)>score(al_neg),
                     score(al), score(al_neg))
    qid <- ifelse(score(al)>score(al_neg),
                  names(unaligned(subject(al))),
                  names(unaligned(subject(al_neg))))
    qid <- do.call(rbind, strsplit(qid, "_"))
    colnames(qid) <- c("dire", "peakID", "tileID")
    ## get max value for qsocre by peakID
    qid <- as.data.frame(qid)
    qid <- cbind(qid, score=qscore)
    qid <- split(qid, qid$peakID)
    qid <- lapply(qid, function(.ele){
      .ele[which.max(.ele$score)[1], ]
    })
    qid <- do.call(rbind, qid)
    .peaks <- .peaks[match(paste(qid$peakID, qid$tileID, sep="_"), .peaks$id)]
    strand(.peaks) <- ifelse(qid$dire=="fwd",
                             ifelse(strand(.peaks)=="-", "-", "+"),
                             ifelse(strand(.peaks)=="-", "+", "-"))
    .peaks$score <- qid$score
    .peaks
  }, subject=subject, query=query, BPPARAM = bpparam)

  peaks <- unlist(GRangesList(peaks))
  ## get Z-score
  m <- mean(peaks$score)
  sd <- sd(peaks$score)
  peaks$Z <- (peaks$score - m)/sd
  S <- length(peaks)
  peaks$pval <- vapply(peaks$score, FUN=function(.ele){
    sum(peaks$score>.ele)/S
  }, FUN.VALUE = numeric(1))
  peaks$adjp <- p.adjust(peaks$pval, method = "BH")

  enhancers(genome = subject@genome, peaks = peaks)
}

#' Output
#' @description Do pairwise alignment for query enhancer to target genome
#' @param query An object of DNAStringSet to represent enhancer
#' @param subject An list of objects of \link{enhancers}.
#' @param \dots Parameters to be used by \link[msa:msa]{msa}.
#' @return An object of \link{enhancers}.
#' @importFrom Biostrings pairwiseAlignment reverseComplement unaligned pattern
#' @importFrom BiocGenerics score
#' @importFrom IRanges subject
#' @importFrom utils combn
#' @importFrom msa msa msaClustalOmega msaClustalW msaMuscle
#' @export
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' hbegfEN <- GRanges("chr14", IRanges(6760805,	6761115))
#' seqEN <- getSeq(BSgenome.Drerio.UCSC.danRer10, hbegfEN)
#' aln_hs <- readRDS(system.file("extdata", "aln_hs.rds",
#'                package="enhancerHomologSearch"))
#' aln_hs$genome <- Hsapiens
#' aln_mm <- readRDS(system.file("extdata", "aln_mm.rds",
#'                package="enhancerHomologSearch"))
#' aln_mm$genome <- Mmusculus
#' al <- alignment(seqEN, list(human=aln_hs, mouse=aln_mm),
#'                 method="ClustalOmega", order="input")
alignment <- function(query, subject, ...){
  stopifnot("query must be an object of DNAStringSet" =
              is(query, "DNAStringSet"))
  stopifnot("The length of query must be 1" = length(query)==1)
  stopifnot("subject must be an list of object of enhancers" =
              is(subject, "list"))
  stopifnot("The length of subject must be more than 0" = length(subject)>0)
  stopifnot("The length of subject must be less than 3" = length(subject)<3)
  null <- lapply(subject, FUN = function(.ele){
    stopifnot("subject must be an list of object of enhancers" =
                is(.ele, "enhancers"))
  })
  if(length(subject)==1){
    q <- getSeq(subject[[1]])
    q <- q[grepl("fwd", names(q))]
    names(q) <- paste(names(subject), names(q), sep = "_")
    if(length(names(query))<1) names(query) <- "Enhancer"
    o <- lapply(seq_along(q), FUN = function(.ele){
      msa(c(query, q[.ele]), ...)
    })
  }else{
    cmb <- combn(x=length(subject), m=2, simplify = FALSE)
    aln <- lapply(cmb, function(x){
      q <- getSeq(subject[[x[1]]])
      s <- getSeq(subject[[x[2]]])
      al <- pairwiseAlignment(rep(q, length(s)),
                              subject = rep(s, each=length(q)),
                              type = "global")
      s1 <- score(al)
      al <- al[order(s1, decreasing = TRUE)]
      s1 <- score(al)
      q <- names(unaligned(pattern(al)))
      s <- names(unaligned(subject(al)))
      id2table <- function(x){
        qid <- do.call(rbind, strsplit(x, "_"))
        colnames(qid) <- c("dire", "ID")
        ## get max value for qsocre by peakID
        qid <- as.data.frame(qid)
      }
      q <- id2table(q)
      s <- id2table(s)
      ids <- cbind(q, s, score=s1)
      al <- al[!duplicated(paste(ids[, 2], ids[, 4]))]
      al
    })
    names(aln) <- vapply(cmb, FUN = function(x){
      paste(names(subject)[x], collapse = "___")
    }, FUN.VALUE = character(1))
    if(length(aln)==1){
      q <- unaligned(pattern(aln[[1]]))
      s <- unaligned(subject(aln[[1]]))
      names(q) <- sub("fwd|rev", sub("___.*$", "", names(aln)), names(q))
      names(s) <- sub("fwd|rev", sub("^.*___", "", names(aln)), names(s))
      if(length(names(query))<1) names(query) <- "Enhancer"
      o <- lapply(seq_along(q), FUN = function(.ele){
        msa(c(query, q[.ele], s[.ele]), ...)
      })
    }else{
      # not support yet.
      stop("More than two species is not supported yet.")
    }
  }
}
