#' Pre-clustered motifs from human and mouse
#'
#' The data were extracted from MotifDb package (v 1.34.0) and grouped by
#' motifStack package (v 1.37.2). The data were packaged as PFMatrixList
#' object by TFBSTools (v 1.30.0)
#'
#' @name motifs
#' @docType data
#' @format a list of PFMatrixList. The names of the list is the grouop distance.
#' @source MotifDb package. Source code for the data generation is in
#' extdata folder
#' @usage data(motifs)
#' @keywords datasets
#' @examples
#'
#' data(motifs)
#' names(motifs)
#' motifs[[1]]
"motifs"
