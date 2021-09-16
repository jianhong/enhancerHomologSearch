#' @import methods
#' @importClassesFrom BSgenome BSgenome
#' @importClassesFrom S4Vectors Annotated
setClassUnion("BSgenomeOrNULL", c("BSgenome", "NULL"))
#' Class \code{"Enhancers"}
#' @description An object of class "Enhancers"
#' represents the output of function \link{getENCODEdata},
#' which includes the sequences of enhancers and their genomic coordinates.
#' @aliases Enhancers
#' @rdname Enhancers-class
#' @slot genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @slot peaks An object of \link[GenomicRanges:GRanges-class]{GRanges}.
#' @slot TFBP An object of \link[Matrix:lsparseMatrix-classes]{lgCMatrix}.
#' @slot TFBP0 An vector of logical.
#' @import methods
#' @importFrom Biostrings DNAStringSet
#' @export
#' @examples
#' Enhancers()
#'
setClass("Enhancers",
         representation(genome = "BSgenomeOrNULL",
                        peaks = "GRanges",
                        TFBP = "lgCMatrix",
                        TFBP0 = "logical"),
         prototype(genome = NULL,
                   peaks = GRanges()),
         validity = function(object){
           if(length(object@peaks)>0 && length(object@genome)>0){
             if(length(intersect(seqlevels(object@peaks),
                                 seqlevels(object@genome)))<1){
               return(
                 "Enhancer seqlevels are not identical.")
             }
           }
           if(length(object@TFBP0) != ncol(object@TFBP)){
             return("TFBPS are not correct.")
           }
           return(TRUE)
         })

#' @rdname Enhancers-class
#' @param genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @param peaks An object of \link[GenomicRanges:GRanges-class]{GRanges}.
#' @param TFBP An object of \link[Matrix:lsparseMatrix-classes]{lgCMatrix}.
#' @param TFBP0 An vector of logical.
#' \code{"Enhancers"}
#' @importFrom Matrix Matrix
#' @export
#' @return An object of Enhancers.
Enhancers <- function(genome, peaks, TFBP, TFBP0){
  if(missing(genome)) genome <- NULL
  if(missing(peaks)) peaks <- GRanges()
  if(missing(TFBP)) TFBP <- new("lgCMatrix")
  if(missing(TFBP0)) TFBP0 <- logical()
  new("Enhancers", genome=genome, peaks=peaks, TFBP=TFBP, TFBP0=TFBP0)
}

#' @rdname Enhancers-class
#' @aliases $,Enhancers-method
#' @param x An object of \code{"Enhancers"}
#' @param name Slot name.
#' @exportMethod `$`
setMethod("$", "Enhancers", function(x, name) slot(x, name))

#' @rdname Enhancers-class
#' @aliases $<-,Enhancers-method
#' @param value The values.
#' @exportMethod `$<-`
setReplaceMethod("$",
                 signature(x="Enhancers"),
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })

#' @rdname Enhancers-class
#' @aliases distance
#' @aliases distance,Enhancers-method
#' @exportMethod `distance`
setMethod("distance", "Enhancers", function(x) slot(x, "peaks")$distance)


if(!exists("distance<-")){
  setGeneric("distance<-", function(x, value) standardGeneric("distance<-"))
}
#' @rdname Enhancers-class
#' @aliases distance<-
#' @aliases distance<-,Enhancers-method
#' @aliases distance<-,Enhancers,ANY-method
#' @exportMethod `distance<-`
setReplaceMethod("distance",
                 signature(x="Enhancers"),
                 function(x, value){
                   x@peaks$distance <- value
                   x
                 })

if(!exists("tfbp")){
  setGeneric("tfbp", function(x) standardGeneric("tfbp"))
}
#' @rdname Enhancers-class
#' @aliases tfbp
#' @aliases tfbp,Enhancers-method
#' @aliases tfbp,Enhancers,ANY-method
#' @exportMethod `tfbp`
setMethod("tfbp",
          signature(x="Enhancers"),
          function(x){
            slot(x, "TFBP")
          })

if(!exists("query_tfbp")){
  setGeneric("query_tfbp", function(x) standardGeneric("query_tfbp"))
}
#' @rdname Enhancers-class
#' @aliases query_tfbp
#' @aliases query_tfbp,Enhancers-method
#' @aliases query_tfbp,Enhancers,ANY-method
#' @exportMethod `query_tfbp`
setMethod("query_tfbp",
          signature(x="Enhancers"),
          function(x){
            slot(x, "TFBP0")
          })

#' @rdname Enhancers-class
#' @aliases getSeq,Enhancers-method
#' @importMethodsFrom Biostrings getSeq reverseComplement
#' @export
#' @param \dots parameters can be passed to upstream functions.
setMethod("getSeq",
          signature(x="Enhancers"),
          function(x, ...){
            genome <- x@genome
            dot <- list(...)
            if(length(dot)>0){
              seq <- getSeq(x = genome, ...)
              if(length(dot[[1]]$id)==0){
                dot[[1]]$id <- paste0(seqnames(dot[[1]]),":",
                                      start(dot[[1]]), "-",
                                      end(dot[[1]]), ":",
                                      strand(dot[[1]]))
              }
            }else{
              seq <- getSeq(x = genome, x@peaks)
              dot <- list(x@peaks)
              dot[[1]]$id <- paste0(seqnames(dot[[1]]),":",
                                    start(dot[[1]]), "-",
                                    end(dot[[1]]), ":",
                                    strand(dot[[1]]))
            }
            names(seq) <- paste("fwd", dot[[1]]$id, sep="_")
            seq2 <- reverseComplement(seq)
            names(seq2) <- sub("fwd", "rev", names(seq))
            c(seq, seq2)
          })

subsetTFBP <- function(x){
  keep <- x@peaks$id
  if(all(keep %in% rownames(x@TFBP))){
    x@TFBP <- x@TFBP[keep, , drop=FALSE]
  }
  x
}

#' @rdname Enhancers-class
#' @aliases subsetByOverlpas,Enhancers-method
#' @param ranges,maxgap,minoverlap,type,invert parameters used by
#' \link[GenomicRanges:findOverlaps-methods]{subsetByOverlaps}
#' @importMethodsFrom IRanges subsetByOverlaps
#' @export
setMethod("subsetByOverlaps",
          signature(x="Enhancers"),
          function(x, ranges, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...){
            x@peaks <- subsetByOverlaps(x@peaks, ranges = ranges,
                                        maxgap = maxgap,
                                        minoverlap = minoverlap,
                                        type = type,
                                        invert = invert,
                                        ...)
            x <- subsetTFBP(x)
            x
          })
#' @rdname Enhancers-class
#' @aliases subset,Enhancers-method
#' @export
setMethod("subset",
          signature(x="Enhancers"),
          function(x, ...){
            x@peaks <- subset(x@peaks, ...)
            x <- subsetTFBP(x)
            x
          })

if(!exists("genome")){
  setGeneric("genome", function(x) standardGeneric("genome"))
}
if(!exists("genome<-")){
  setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))
}

#' @rdname Enhancers-class
#' @aliases seqinfo,Enhancers-method
#' @importFrom GenomeInfoDb seqinfo
#' @importMethodsFrom GenomeInfoDb seqinfo
#' @exportMethod `seqinfo`
setMethod("seqinfo",
          signature(x="Enhancers"),
          function(x){
            seqinfo(x@genome)
          })
#' @rdname Enhancers-class
#' @aliases genome
#' @aliases genome,Enhancers-method
#' @importFrom GenomeInfoDb genome
#' @importMethodsFrom GenomeInfoDb genome
#' @exportMethod `genome`
setMethod("genome",
         signature(x="Enhancers"),
         function(x){
           x@genome
         })
#' @rdname Enhancers-class
#' @aliases genome<-
#' @aliases genome<-,Enhancers-method
#' @aliases genome<-,Enhancers,BSgenome-method
#' @importFrom GenomeInfoDb genome<-
#' @importMethodsFrom GenomeInfoDb genome<-
#' @exportMethod `genome<-`
setReplaceMethod("genome",
                 signature(x="Enhancers"),
                 function(x, value){
                   x@genome <- value
                   x
                 })


if(!exists("peaks")){
  setGeneric("peaks", function(x) standardGeneric("peaks"))
}
if(!exists("peaks<-")){
  setGeneric("peaks<-", function(x, value) standardGeneric("peaks<-"))
}

#' @rdname Enhancers-class
#' @aliases peaks
#' @aliases peaks,Enhancers-method
#' @exportMethod `peaks`
setMethod("peaks",
          signature(x="Enhancers"),
          function(x){
            x@peaks
          })
#' @rdname Enhancers-class
#' @aliases peaks<-
#' @aliases peaks<-,Enhancers-method
#' @aliases peaks<-,Enhancers,GRanges-method
#' @exportMethod `peaks<-`
setReplaceMethod("peaks",
                 signature(x="Enhancers"),
                 function(x, value){
                   x@peaks <- value
                   x
                 })

#' @name coerce
#' @rdname Enhancers-class
#' @aliases coerce,Enhancers,GRanges-method
#' @exportMethod coerce
setAs(from="Enhancers", to="GRanges", function(from){
  from@peaks
})
#' @rdname Enhancers-class
#' @aliases show,Enhancers-method
#' @param object An object of \code{"Enhancers"}
#' @export
setMethod("show",
          signature(object="Enhancers"),
          function(object){
            if(length(object@peaks)){
              cat("This is an object with ",
                  length(object@peaks),
                  " Enhancers for ",
                  organism(object@genome))
            }else{
              cat("This is an empty object of Enhancers.")
            }
          })
