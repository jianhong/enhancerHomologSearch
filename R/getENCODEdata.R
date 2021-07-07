#' Download enhancer sequences from ENCODE
#' @description Query enhancer peak and extract sequences from ENCODE
#' @param genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @param markers Enhancer markers. Default 'H3K4me1'. For active enhancer,
#' it can be set as c('H3K4me1', 'H3K27ac')
#' @param window_size,step The size of windows and steps to split the peaks
#' into small pieces. These parameter is used because the width of
#' histone marker peaks are different sizes. Break the peaks into small
#' pieces can increase the matching score and align the matching for
#' different peaks into same size. The window_size is also be used for
#' overlapping detection of multiple histone markers.
#' @param \dots Parameters can be passed to \link{queryEncode}
#' @return An object of \link{enhancers} with genome, and peaks.
#' The peaks is an object of GRanges. The genome is an object of BSgenome.
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
#'                     partialMatch=c(biosample_summary="spinal cord"))
getENCODEdata <- function(genome,
                          markers="H3K4me1",
                          window_size = 1000L,
                          step = 50L,
                          ...) {
  stopifnot('genome must be an object of BSgenome'=is(genome, "BSgenome"))
  dir = tempdir() #The directory for saving peak files.
  org <- organism(genome)
  assembly <- guessAssembly(genome)
  names(assembly) <- rep("assembly", length(assembly))
  names(markers) <- rep("target.label", length(markers))
  args <- list(...)
  if(!"exactMatch" %in% names(args)){
    exactMatch <- markers
  }else{
    exactMatch <- exactMatch[!names(exactMatch) %in% "target.label"]
    exactMatch <- c(markers, exactMatch)
  }
  if(!"assembly" %in% names(exactMatch)){
    exactMatch <- c(exactMatch, assembly)
  }
  if(!"replicates.library.biosample.donor.organism.scientific_name" %in%
     names(exactMatch)){
    exactMatch <-
      c(replicates.library.biosample.donor.organism.scientific_name=org,
        exactMatch)
  }
  res <- queryEncode(exactMatch=exactMatch, ...)

  if(length(res)==0){
    stop("No data available by current setting.")
  }

  datasets <- vapply(res, FUN = function(.ele) .ele$`@id`,
                     FUN.VALUE = character(1))
  names(datasets) <- rep("dataset", length(datasets))

  args$exactMatch <- c(type='File', output_type="peaks",
                       file_format="bed", assembly,
                       datasets)

  res2 <- do.call(queryEncode, args = args)
  if(length(res2)==0){
    stop("No data available by current setting.")
  }

  urls <- lapply(res2, FUN = function(.ele) .ele$cloud_metadata$url)
  format <- lapply(res2, FUN = function(.ele) .ele$file_format_type)

  peaks <- mapply(import, urls, format, SIMPLIFY = FALSE)
  peaks <- split(peaks,
                 vapply(res2, `[[`, FUN.VALUE = character(1), i="target"))
  peaks <- lapply(peaks, function(.ele) reduce(unlist(GRangesList(.ele))))
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
