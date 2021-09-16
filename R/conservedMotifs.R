#' check the conserved motifs in the orthologs
#'
#' Print the conserved motifs in the alignments
#'
#' @param aln alignment of multiple DNAs. Output of \link{alignment} function.
#' @param aln_hs,aln_mm output of link{searchTFBPS} for human and mouse.
#' @param PWMs The Position Weight Matrix list represented as a numeric matrix.
#' Object of \link[TFBSTools:XMatrixList]{PWMatrixList} or
#' \link[TFBSTools:XMatrixList]{PFMatrixList}.
#' @param queryGenome An object of \link[BSgenome:BSgenome-class]{BSgenome} for
#' query enhancer.
#' @param background background nucleotide frequencies. Default is "genome".
#' Refer \link[motifmatchr]{matchMotifs} for details.
#' @importFrom motifmatchr matchMotifs
#' @importFrom Biostrings unmasked mask
#' @importFrom IRanges gaps
#' @importMethodsFrom Matrix t rowSums
#' @export
#' @return A list of XStringviews.
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
#' data(motifs)
#' conservedMotifs(al[[1]], aln_hs, aln_mm, motifs[["dist60"]], Drerio)
conservedMotifs <- function(aln, aln_hs, aln_mm, PWMs,
                            queryGenome, background="genome"){
  stopifnot(is(aln, "DNAMultipleAlignment"))
  stopifnot(is(queryGenome, "BSgenome"))
  if(!missing(aln_hs)) stopifnot(is(aln_hs, "Enhancers"))
  if(!missing(aln_mm)) stopifnot(is(aln_mm, "Enhancers"))
  background <- match.arg(background, choices = c("subject", "genome", "even"))
  TFBPS_target <- 1
  ## get the matched motifs for query enhancer
  if(!missing(aln_hs)){
    TFBPS_target <- TFBPS_target & query_tfbp(aln_hs)
  }
  if(!missing(aln_mm)){
    TFBPS_target <- TFBPS_target & query_tfbp(aln_mm)
  }
  seq <- unmasked(aln)
  n <- names(seq)
  TFBPS <- TFBPS_target
  genome <- queryGenome
  if(!missing(aln_hs)){
    h <- n[grepl("human", n)]
    TFBPS_human <- tfbp(subsetByOverlaps(aln_hs,
                                         GRanges(sub("human_", "", h)),
                                         type = "equal"))
    if(length(dim(TFBPS_human))){
      if(nrow(TFBPS_human)==1){
        TFBPS <- rbind(TFBPS, TFBPS_human)
        genome <- c(genome, genome(aln_hs))
      }
    }
  }
  if(!missing(aln_mm)){
    m <- n[grepl("mouse", n)]
    TFBPS_mouse <- tfbp(subsetByOverlaps(aln_mm,
                                         GRanges(sub("mouse_", "", m)),
                                         type = "equal"))
    if(length(dim(TFBPS_mouse))){
      if(nrow(TFBPS_mouse)==1){
        TFBPS <- rbind(TFBPS, TFBPS_mouse)
        genome <- c(genome, genome(aln_mm))
      }
    }
  }
  if(length(dim(TFBPS))!=2){
    stop("no data available")
  }
  TFBPS <- t(TFBPS)
  PWMs <- PWMs[rownames(TFBPS[rowSums(TFBPS)==ncol(TFBPS), , drop=FALSE])]
  if(length(PWMs)<1){
    stop("no data available")
  }
  seq <- gsub("-", "", seq)
  mm <- mapply(seq, genome, FUN=function(s, g){
    matchMotifs(PWMs, subject = s, out="positions", genome=g, bg=background)
  }, SIMPLIFY = FALSE)
  mm <- swapList(mm)
  mm <- lapply(mm, function(.ele){
    lapply(.ele, unlist)
  })
  mm <- mm[vapply(mm, FUN=function(.ele) all(lengths(.ele)>0),
                  FUN.VALUE = logical(1))]
  if(length(mm)==0){
    stop("no data available")
  }
  mm <- swapList(mm)
  mm <- lapply(mm, function(.ele) lapply(.ele, reduce))
  mm <- lapply(mm, IRangesList)
  mm <- lapply(mm, unlist)
  seq <- DNAStringSet(seq)
  gaps <- mapply(gaps, mm, 1, lengths(seq))
  mapply(function(.seq, .gap) mask(.seq, start(.gap), end(.gap)),
         seq, gaps)
}
