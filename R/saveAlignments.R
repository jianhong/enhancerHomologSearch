checkMotifCons <- function(motifConsensus){
  if(length(motifConsensus)){
    stopifnot("motifConsensus must be an object of DNAStringSet"=
                is(motifConsensus, "DNAStringSet"))
    af <- alphabetFrequency(motifConsensus)
    af <- colSums(af)
    stopifnot("motifConsensus can only contain 'A', 'C', 'G', 'T', 'N'"=
                sum(lengths(motifConsensus))==
                sum(af[c("A", "C", "G", "T", "N")]))
  }
}
#' @importFrom Biostrings consensusMatrix DNAStringSet write.phylip
#' matchPattern alphabetFrequency injectHardMask
#' @importFrom IRanges IRangesList Views
saveAln <- function (x, filepath, motifConsensus) {
  stopifnot("input must be an object of MultipleAlignment"=
              inherits(x, "MultipleAlignment"))
  checkMotifCons(motifConsensus)
  ## get consensus
  enhancer <- strsplit(as.character(x@unmasked[1]), "")[[1]]
  cons <- consensusMatrix(x, baseOnly=TRUE)
  consid <- apply(cons[rownames(cons) %in% c("A", "C", "G", "T"), ,
                       drop=FALSE], 2, which.max)
  consid[apply(cons[rownames(cons) %in% c("A", "C", "G", "T"), ,
                    drop=FALSE], 2, max)<2] <-
    which(!rownames(cons) %in% c("A", "C", "G", "T"))[1]
  cons <- rownames(cons)[consid]
  cons[cons!=enhancer] <- "-"
  cons <- paste(cons, collapse = "")
  cons_m <- gsub("-", "N", cons)
  cons <- DNAStringSet(cons)
  names(cons) <- "Consensus"
  x@unmasked <- c(x@unmasked, cons)
  ## search binding sites
  if(length(motifConsensus)){
    mt <- lapply(motifConsensus, function(.ele)
      matchPattern(pattern = .ele, subject = cons_m))
    mt <- mt[lengths(mt)>0]
    mt <- unlist(IRangesList(mt))
    mt <- reduce(mt)
    if(length(mt)){
      m <- Views(cons_m, start = start(mt), end = end(mt))
      m <- injectHardMask(m, letter="-")
      m <- DNAStringSet(m)
      names(m) <- "motifConsensus"
      x@unmasked <- c(x@unmasked, m)
    }
  }
  write.phylip(x, filepath = filepath)
}
#' output alignments
#' @description Save enhancer homologs to file in phylip format.
#' @param al output of \link{alignment}.
#' @param output_folder output folder.
#' @param motifConsensus Transcription factor binding consensus.
#' @export
#' @return The I/O status.
#' @examples
#' al <- readRDS(system.file("extdata", "al.rds",
#'                package="enhancerHomologSearch"))
#' tmpfolder <- tempdir()
#' library(MotifDb)
#' motifs <- query(MotifDb, "JASPAR_CORE")
#' consensus <- sapply(motifs, consensusString)
#' consensus <- DNAStringSet(gsub("\\?", "N", consensus))
#' saveAlignments(al, output_folder=tmpfolder, motifConsensus=consensus)
saveAlignments <- function(al, output_folder=tempdir(), motifConsensus=NULL){
  checkMotifCons(motifConsensus)
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  if(length(names(al))!=length(al)){
    names(al) <- paste0("aln", seq_along(al), ".phylip.txt")
  }else{
    names(al) <- paste0(names(al), ".phylip.txt")
  }
  null <- mapply(saveAln, x=al,
                 filepath = file.path(output_folder, names(al)),
                 MoreArgs=list(motifConsensus = motifConsensus))
  return(invisible(null))
}
