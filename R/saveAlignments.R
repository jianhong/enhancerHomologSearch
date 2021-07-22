#' @importFrom Biostrings consensusMatrix DNAStringSet write.phylip
saveAln <- function (x, filepath) {
  stopifnot("input must be an object of MultipleAlignment"=
              inherits(x, "MultipleAlignment"))
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
  cons <- DNAStringSet(cons)
  names(cons) <- "Consensus"
  x@unmasked <- c(x@unmasked, cons)
  write.phylip(x, filepath = filepath)
}
#' output alignments
#' @description Save enhancer homologs to file in phylip format.
#' @param al output of \link{alignment}.
#' @param output_folder output folder.
#' @export
#' @return The I/O status.
#' @examples
#' al <- readRDS(system.file("extdata", "al.rds",
#'                package="enhancerHomologSearch"))
#' tmpfolder <- tempdir()
#' saveAlignments(al, output_folder=tmpfolder)
saveAlignments <- function(al, output_folder=tempdir()){
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  if(length(names(al))!=length(al)){
    names(al) <- paste0("aln", seq_along(al), ".phylip.txt")
  }else{
    names(al) <- paste0(names(al), ".phylip.txt")
  }
  null <- mapply(saveAln, x=al, filepath = file.path(output_folder, names(al)))
  return(invisible(null))
}
