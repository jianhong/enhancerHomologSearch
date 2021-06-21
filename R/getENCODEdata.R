#' Download enhancer sequences from ENCODE
#' @description Query enhancer peak and extract sequences from ENCODE
#' @param genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @param markers Enhancer markers. Default 'H3K4me1'. For active enhancer,
#' it can be set as c('H3K4me1', 'H3K27ac')
#' @param file_format File format for the peak files. Default is bed.
#' @param output_filter The filter for output type of peak files.
#' Default is "stable peaks".
#' @param window_size,step The size of windows and steps to split the peaks
#' into small pieces. These parameter is used because the width of
#' histone marker peaks are different sizes. Break the peaks into small
#' pieces can increase the matching score and align the matching for
#' different peaks into same size. The window_size is also be used for
#' overlapping detection of multiple histone markers.
#' @param \dots The parameters could be used by
#' \link[ENCODExplorer:queryEncode]{queryEncode}
#' @return An object of \link{enhancers} with genome, and peaks.
#' The peaks is an object of GRanges. The genome is an object of BSgenome.
#' @importFrom ENCODExplorer queryEncode downloadEncode
#' @importFrom rtracklayer import
#' @importFrom Biostrings getSeq
#' @importFrom BiocGenerics organism
#' @importFrom IRanges subsetByOverlaps
#' @import methods
#' @import GenomicRanges
#' @export
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' hs <- getENCODEdata(genome=Hsapiens,
#'                     biosample_name = "spinal cord",
#'                     biosample_type = "tissue")
getENCODEdata <- function(genome,
                          markers = c("H3K4me1"),
                          file_format = "bed",
                          output_filter = c("stable peaks", "replicated peaks"),
                          window_size = 1000L,
                          step = 50L,
                          ...) {
  stopifnot('genome must be an object of BSgenome'=is(genome, "BSgenome"))
  dir = tempdir() #The directory for saving peak files.
  org <- organism(genome)
  assembly <- guessAssembly(genome)
  peaks <- lapply(markers, FUN=function(marker){
    q <- queryEncode(organism = org,
                     target = marker,
                     file_format = file_format,
                     ...)
    keep <- q$output_type %in% output_filter &
      q$assembly %in% assembly
    if(sum(keep)<1){
      message("The output of searching results before filtering:")
      print(q)
      stop("No data available by current setting:",
           "output_filter == ", output_filter,
           " and assembly == ", paste(assembly, collapse="/"),
           "Available output_filter: ", unique(q$output_type),
           "Avallable assembly: ", unique(q$assembly))
    }
    q <- q[keep, , drop = FALSE]
    f <- downloadEncode(q, dir = dir)
    p <- tryCatch({
        if(file_format=="bed"){
          lapply(f, import, format = "narrowPeak")
        }else{
          lapply(f, import)
        }
      },
      error = function(e){
        warning(e)
      })
    unlink(f)
    p <- reduce(unlist(GRangesList(p)))
  })
  ## get intersection of the ranges
  subsetByOverlapsWindow <- function(a, b){
    subsetByOverlaps(a, b, maxgap = window_size)
  }
  if(length(peaks)>1){
    peaks <- Reduce(f = subsetByOverlapsWindow, x = peaks)
  }else{
    peaks <- peaks[[1]]
  }
  if(length(peaks)<2){
    stop("There is less than 2 peaks left. Please try to reduce the markers")
  }
  ## reset the peak ranges
  w <- width(peaks)
  w1 <- ceiling(ceiling(w/window_size) * window_size / 2)
  start(peaks) <- end(peaks) <- ceiling((start(peaks) + end(peaks))/2)
  start(peaks) <- start(peaks) - w1
  width(peaks) <- 2 * w1
  peaks <- slidingWindows(peaks, width = window_size, step = step)
  pid <- rep(seq_along(peaks), lengths(peaks))
  pid2 <- unlist(lapply(lengths(peaks), FUN = seq.int))
  peaks <- unlist(peaks)
  peaks$id <- paste(pid, pid2, sep="_")
  #peaks <- peaks[width(peaks)==window_size]
  enhancers(genome = genome, peaks = peaks)
}
