#' shuffle reads
#'
#' Uses the uShuffle library to shuffle reads
#'
#' @param reads An object of \link[Biostrings:XStringSet-class]{BStringSet}.
#' @param k the k-let size.
#' @param n the number of random sequences to generate.
#' @return An object of \link[Biostrings:XStringSet-class]{BStringSet}.
#' @importFrom Biostrings readDNAStringSet DNAStringSet RNAStringSet
#'  AAStringSet
#' @importFrom methods is
#' @export
#' @references Jiang, M., Anderson, J., Gillespie, J. et al.
#' uShuffle: A useful tool for shuffling biological sequences while preserving
#'  the k-let counts. BMC Bioinformatics 9, 192 (2008).
#'  https://doi.org/10.1186/1471-2105-9-192
#' @examples
#' library(Biostrings)
#' f <- DNAStringSet(c("CTC-NACCAGTAT", "TTGA", "TACCTAGAG"))
#' shuffle(f)

shuffle <- function(reads, k=2, n=2){
  if(is.character(reads)){
    reads <- readDNAStringSet(reads)
  }
  stopifnot("reads must be an object of BStringSet"=
              inherits(reads, c("BStringSet", "DNAStringSet", "RNAStringSet",
                                "AAStringSet", "DNAString", "RNAString",
                                "AAString")))
  if(is(reads, "DNAString")){
    reads <- DNAStringSet(reads)
  }
  if(is(reads, "RNAString")){
    reads <- RNAStringSet(reads)
  }
  if(is(reads, "AAString")){
    reads <- AAStringSet(reads)
  }
  if(length(reads)==0){
    return(NULL)
  }
  in_seqs <- as.character(reads)
  out_seqs <- rushuffle(in_seqs, k, n)
  if(is(reads, "DNAStringSet")){
    seqs <- DNAStringSet(out_seqs)
  }
  if(is(reads, "RNAStringSet")){
    seqs <- RNAStringSet(out_seqs)
  }
  if(is(reads, "AAStringSet")){
    seqs <- AAStringSet(out_seqs)
  }
  if(length(names(reads))){
    names(seqs) <- paste0(rep(names(reads), each=n),
                          "_shuffle_",
                          rep(seq.int(n), length(reads)))
  }
  return(seqs)
}
