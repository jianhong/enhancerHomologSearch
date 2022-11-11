#' check the conserved motifs in the orthologs
#'
#' Print the conserved motifs in the alignments
#'
#' @param aln alignment of multiple DNAs. Output of \link{alignment} function.
#' @param aln_list The list of output of \link{searchTFBPS}
#'  such as for human and mouse.
#' @param PWMs The Position Weight Matrix list represented as a numeric matrix.
#' Object of \link[TFBSTools:XMatrixList]{PWMatrixList} or
#' \link[TFBSTools:XMatrixList]{PFMatrixList}.
#' @param queryGenome An object of \link[BSgenome:BSgenome-class]{BSgenome} for
#' query enhancer.
#' @param background Background nucleotide frequencies. Default is "genome".
#' Refer \link[motifmatchr]{matchMotifs} for details.
#' @param \dots Other parameters can be passed to to
#' \link[motifmatchr:matchMotifs]{matchMotifs}.
#' @param output_folder Output folder name.
#' @param format The format of output files with motif match positions.
#'  Available formats are 'txt' and 'html'. Default is 'txt'.
#' @importFrom motifmatchr matchMotifs
#' @importFrom Biostrings unmasked mask getSeq DNAMultipleAlignment
#' reverseComplement
#' @importFrom IRanges gaps Views
#' @importMethodsFrom Matrix t rowSums
#' @export
#' @return A list of \link[Biostrings:XStringViews-class]{XStringViews}
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
#' conservedMotifs(al[[1]], list(human=aln_hs, mouse=aln_mm),
#'                 motifs[["dist60"]], Drerio)
conservedMotifs <- function(aln, aln_list, PWMs,
                            queryGenome,
                            background="genome",
                            ...,
                            output_folder,
                            format=c("txt", "html")){
  stopifnot(is(aln, "DNAMultipleAlignment"))
  stopifnot(is(queryGenome, "BSgenome"))
  checkSubject(aln_list, subjectIsList=TRUE)
  background <- match.arg(background, choices = c("subject", "genome", "even"))
  ## get sequences of aligned enhancers
  seq <- unmasked(aln)
  n <- names(seq)[-1]
  n2 <- names(aln_list)
  n <- sub("_.*$", "", n)
  stopifnot("The aln_list should keep same as input for alignment"=
              all(n %in% n2) && all(n2 %in% n))
  seq <- c(seq[1], seq[-1][match(n2, n)])
  n <- names(seq)
  genome <- queryGenome
  ### get hiting motifs
  TFBPS <- mapply(aln_list, n2, FUN=function(enh, prefix){
    h <- n[grepl(prefix, n)]
    tfbps <- tfbp(subsetByOverlaps(enh,
                                   GRanges(sub(paste0(prefix, "_"), "", h)),
                                   type = "equal"))
    if(length(dim(tfbps))){
      if(nrow(tfbps)==1){
        return(tfbps)
      }
    }
    stop("Something not supported. Please try to run `searchTFBPS`.")
  })
  TFBPS <- do.call(rbind, TFBPS)
  genome <- lapply(aln_list, genome)
  genome <- c(queryGenome, genome)

  if(length(dim(TFBPS))!=2){
    stop("no data available: TFBPS is not a two dim data.")
  }

  TFBPS <- t(TFBPS)
  pwms <- rownames(TFBPS[rowSums(TFBPS)==ncol(TFBPS), , drop=FALSE])
  PWMs <- PWMs[names(PWMs) %in% pwms]
  if(length(PWMs)<1){
    stop("no data available: not cosvered motifs in all inputs.")
  }
  PWMs_rev <- lapply(PWMs, reverseComplement)
  names(PWMs_rev) <- paste0("fwd_", names(PWMs))
  PWMs_rev <- c(PWMs, PWMs_rev)
  seq <- gsub("-", "", seq)
  mm <- mapply(seq, genome, FUN=function(s, g){
    matchMotifs(PWMs_rev, subject = s, out="positions",
                genome=g, bg=background,
                ...)
  }, SIMPLIFY = FALSE)
  mm <- swapList(mm)
  mm <- lapply(mm, function(.ele){
    lapply(.ele, unlist)
  })
  mm <- mm[vapply(mm, FUN=function(.ele) all(lengths(.ele)>0),
                  FUN.VALUE = logical(1))]
  if(length(mm)==0){
    stop("no data available: no matched consensus.")
  }
  mm <- swapList(mm)
  mm <- lapply(mm, function(.ele){
    names(.ele) <- sub("fwd_", "", names(.ele))
    .ele <- split(.ele, names(.ele))
    .ele <- lapply(.ele, IRangesList)
    .ele <- lapply(.ele, unlist)
  })
  mm <- lapply(mm, function(.ele) lapply(.ele, reduce))
  mm <- lapply(mm, IRangesList)
  mm <- lapply(mm, unlist)
  seq <- DNAStringSet(seq)
  view_list <- mapply(Views, seq, mm)
  if(!missing(output_folder)){
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    involved_motifs <- as(IRangesList(mm), "GRanges")
    involved_motifs$motif <- unlist(lapply(mm, names))
    involved_motifs$seq <- getSeq(seq, involved_motifs)
    write.csv(involved_motifs,
              file = file.path(output_folder, "motifhits.csv"),
              row.names=FALSE)
    mapply(seq, mm, names(seq), FUN=function(.seq, .match, .name){
      v <- lapply(seq_along(.match), function(i)
        Views(subject=.seq, start = .match[i]))
      names(v) <- names(.match)
      v <- lapply(v, injectHardMask, letter="-")
      a <- DNAStringSet(.seq)
      names(a) <- .name
      v <- c(a, DNAStringSet(v))
      x <- DNAMultipleAlignment(v)
      filepath <- file.path(output_folder,
                            paste(make.names(
                              gsub(":", "_",
                                   sub("[+*]$", "p",
                                       sub("-$", "n", .name)))),
                                  format, sep="."))
      if(format=="txt"){
        write.phylip(x, filepath = filepath)
      }else{
        write.phylip_html(x, filepath = filepath)
      }
    })
    return(invisible(view_list))
  }else{
    view_list
  }
}
