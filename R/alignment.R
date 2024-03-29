#' Get alignment scores
#' @description Do pairwise alignment for query enhancer to target genome
#' @param query An object of DNAStringSet to represent enhancer
#' @param subject Output of getENCODEdata. An object of \link{Enhancers}
#' @param block The size of sequences to do alignment. Increase the size will
#' increase the memory cost. Default 1000.
#' @param bpparam BiocParallel parameters.
#' @param \dots not used.
#' @return An object of \link{Enhancers}.
#' @importFrom Biostrings reverseComplement unaligned
#' @importFrom pwalign pairwiseAlignment
#' @importFrom BiocGenerics score
#' @importFrom IRanges subject
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
#' @export
#' @examples
#' library(BiocParallel)
#' bpparam <- MulticoreParam(workers = 1, tasks=200, progressbar=TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' peaks <- GRanges("chr1", IRanges(seq(5000, 50000, by=1000), width=1000))
#' peaks$id <- paste(seq_along(peaks), 1, sep="_")
#' subj <- Enhancers(genome=Hsapiens, peaks=peaks)
#' q <- getSeq(Hsapiens, GRanges("chr1", IRanges(90000, width=1000)))
#' ao <- alignmentOne(q, subj, bpparam=bpparam)
alignmentOne <- function(query, subject, block=1000, bpparam = bpparam(), ...){
  checkQuerySubject(query, subject, subjectIsList=FALSE)
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
    stopifnot("subject must be output of getENCODEdata."=ncol(qid)==3)
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

  Enhancers(genome = subject@genome, peaks = peaks)
}

#' Output
#' @description Do pairwise alignment for query enhancer to target genome
#' @param query An object of DNAStringSet to represent enhancer
#' @param subject An list of objects of \link{Enhancers}.
#' @param method specifies the multiple sequence alignment to be used;
#' currently, "ClustalW", and "Muscle" are supported. Default is "Muscle"
#' @param cluster The clustering method which should be used.
#' Possible values are "nj" (default) and "upgma".
#' In the original ClustalW implementation, this parameter is called clustering.
#' @param gapOpening gap opening penalty; the default is 400 for DNA sequences and 420 for RNA sequences. The default for amino acid sequences depends on the profile score settings: for the setting le=TRUE, the default is 2.9, for sp=TRUE, the default is 1,439, and for sv=TRUE, the default is 300. Note that these defaults may not be suitable if custom substitution matrices are being used. In such a case, a sensible choice of gap penalties that fits well to the substitution matrix must be made.
#' @param gapExtension gap extension penalty; the default is 0.
#' @param maxiters maximum number of iterations; the default is 16.
#' @param substitutionMatrix substitution matrix for scoring matches and
#' mismatches; The valid choices for this parameter are "iub" and "clustalw".
#' In the original ClustalW implementation, this parameter is called matrix.
#' @param order how the sequences should be ordered in the output object; if "aligned" is chosen, the sequences are ordered in the way the multiple sequence alignment algorithm orders them. If "input" is chosen, the sequences in the output object are ordered in the same way as the input sequences.
#' @param \dots Parameters can be used by Muscle, or ClustalW.
#' @return An object of \link{Enhancers}.
#' @importFrom Biostrings reverseComplement unaligned pattern
#' @importFrom pwalign pairwiseAlignment
#' @importFrom BiocGenerics score
#' @importFrom IRanges subject
#' @importFrom utils combn
#' @import Rcpp
#' @export
#' @useDynLib enhancerHomologSearch
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(BSgenome.Drerio.UCSC.danRer10)
#' LEN <- GRanges("chr4", IRanges(19050041, 19051709))
#' seqEN <- getSeq(BSgenome.Drerio.UCSC.danRer10, LEN)
#' aln_hs <- readRDS(system.file("extdata", "aln_hs.rds",
#'                package="enhancerHomologSearch"))
#' genome(aln_hs) <- Hsapiens
#' aln_mm <- readRDS(system.file("extdata", "aln_mm.rds",
#'                package="enhancerHomologSearch"))
#' genome(aln_mm) <- Mmusculus
#' al <- alignment(seqEN, list(human=aln_hs, mouse=aln_mm),
#'                 method="ClustalW", order="input")
alignment <- function(query, subject,
                      method=c("ClustalW", "Muscle"),
                      cluster=c("nj", "upgma", "upgmamax", "upgmamin", "upgmb"),
                      substitutionMatrix=c("iub", "clustalw"),
                      gapOpening=ifelse(method[1]=="ClustalW", 15.0, 400),
                      gapExtension=ifelse(method[1]=="ClustalW", 6.66, 0),
                      maxiters=ifelse(method[1]=="ClustalW", 3, 16),
                      order=c("aligned", "input"),
                      ...){
  checkQuerySubject(query, subject, subjectIsList=TRUE)
  args <- as.list(match.call(expand.dots=FALSE))
  args$query <- NULL
  args$subject <- NULL
  args$cluster <- match.arg(cluster)
  args$method <- match.arg(method)
  args$order <- match.arg(order)
  if(is.character(substitutionMatrix)){
    args$substitutionMatrix <- match.arg(substitutionMatrix)
  }
  args$gapExtension <- eval(quote(gapExtension))
  args$gapOpening <- eval(quote(gapOpening))
  args$maxiters <- eval(quote(maxiters))
  if(length(subject)==1){
    q <- getSeq(subject[[1]])
    q <- q[grepl("fwd", names(q))]
    names(q) <- paste(names(subject), names(q), sep = "_")
    if(length(names(query))<1) names(query) <- "Enhancer"
    o <- lapply(seq_along(q), FUN = function(.ele){
      args$inputSeqs <- c(query, q[.ele])
      do.call(msa, args)
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
    if(length(names(query))<1) names(query) <- "Enhancer"
    if(length(aln)==1){
      q <- unaligned(pattern(aln[[1]]))
      s <- unaligned(subject(aln[[1]]))
      names(q) <- sub("fwd|rev", sub("___.*$", "", names(aln)), names(q))
      names(s) <- sub("fwd|rev", sub("^.*___", "", names(aln)), names(s))
      o <- lapply(seq_along(q), FUN = function(.ele){
        args$inputSeqs <- c(query, q[.ele], s[.ele])
        do.call(msa, args)
      })
    }else{
      aln <- mapply(aln,
                    sub("___.*$", "", names(aln)),
                    sub("^.*___", "", names(aln)),
                    FUN=function(.ele, .qn, .sn){
                      q <- unaligned(pattern(.ele))
                      s <- unaligned(subject(.ele))
                      names(q) <- sub("fwd|rev", .qn, names(q))
                      names(s) <- sub("fwd|rev", .sn, names(s))
                      df <- data.frame(q=names(q), s=names(s))
                      colnames(df) <- c(.qn, .sn)
                      list(df=df, seq=c(q, s))
                    }, SIMPLIFY=FALSE)
      seqs <- lapply(aln, `[[`, i="seq")
      seqs <- Reduce(c, seqs)
      seqs <- unique(seqs)
      aln <- lapply(aln, `[[`, i="df")
      aln <- Reduce(merge, aln)
      o <- apply(aln, 1, FUN = function(.ele){
        args$inputSeqs <- c(query, seqs[.ele])
        do.call(msa, args)
      }, simplify = FALSE)
    }
  }
}
