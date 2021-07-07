#' @import methods
#' @importClassesFrom BSgenome BSgenome
#' @importClassesFrom S4Vectors Annotated
setClassUnion("BSgenomeOrNULL", c("BSgenome", "NULL"))
#' Class \code{"enhancers"}
#' @description An object of class "enhancers"
#' represents the output of function \link{getENCODEdata},
#' which includes the sequences of enhancers and their genomic coordinates.
#' @aliases enhancers
#' @rdname enhancers-class
#' @slot genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @slot peaks An object of \link[GenomicRanges:GRanges-class]{GRanges}.
#' @import methods
#' @importFrom Biostrings DNAStringSet
#' @export
#' @examples
#' enhancers()
#'
setClass("enhancers",
         representation(genome = "BSgenomeOrNULL",
                        peaks = "GRanges"),
         prototype(genome = NULL,
                   peaks = GRanges()),
         validity = function(object){
           if(length(object@peaks)>0 && length(object@genome)>0){
             if(length(intersect(seqlevels(object@peaks),
                                 seqlevels(object@genome)))<1){
               return(
                 "Please use getENCODEdata to create enhancers object.")
             }
           }
           return(TRUE)
         })

#' @rdname enhancers-class
#' @param \dots Each argument in \dots becomes an slot in the new("enhancers")
#' or An object of \code{"GRanges"} for getSeq,enhancers-method
#' \code{"enhancers"}
#' @export
#' @return An object of enhancers.
enhancers <- function(...){
  new("enhancers", ...)
}

#' @rdname enhancers-class
#' @aliases $,enhancer-method
#' @aliases $<-,enhancer-method
#' @param x An object of \code{"enhancers"}
#' @param name Slot name.
#' @exportMethod `$`
setMethod("$", "enhancers", function(x, name) slot(x, name))

#' @rdname enhancers-class
#' @aliases $<-,enhancer-method
#' @param value The values.
#' @exportMethod `$<-`
setReplaceMethod("$",
                 signature(x="enhancers"),
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })
#' @rdname enhancers-class
#' @aliases getSeq,enhancers-method
#' @importMethodsFrom Biostrings getSeq reverseComplement
#' @export
setMethod("getSeq",
          signature(x="enhancers"),
          function(x, ...){
            genome <- x@genome
            dot <- list(...)
            if(length(dot)>0){
              seq <- getSeq(x = genome, ...)
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
#' @rdname enhancers-class
#' @aliases subsetByOverlpas,enhancers-method
#' @param ranges,maxgap,minoverlap,type,invert parameters used by
#' \link[GenomicRanges:findOverlaps-methods]{subsetByOverlaps}
#' @importMethodsFrom IRanges subsetByOverlaps
#' @export
setMethod("subsetByOverlaps",
          signature(x="enhancers"),
          function(x, ranges, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...){
            x@peaks <- subsetByOverlaps(x@peaks, ranges = ranges,
                                        maxgap = maxgap,
                                        minoverlap = minoverlap,
                                        type = type,
                                        invert = invert,
                                        ...)
            x
          })
#' @rdname enhancers-class
#' @aliases subset,enhancers-method
#' @export
setMethod("subset",
          signature(x="enhancers"),
          function(x, ...){
            x@peaks <- subset(x@peaks, ...)
            x
          })
#' @name genome
#' @rdname enhancers-class
#' @aliases genome,enhancers-method
#' @importFrom GenomeInfoDb genome
#' @importMethodsFrom GenomeInfoDb genome
#' @export
if(!exists("genome")){
  setGeneric("genome", function(x) standardGeneric("genome"))
}
setMethod("genome",
         signature(x="enhancers"),
         function(x){
           x@genome
         })
#' @name coerce
#' @rdname enhancers-class
#' @aliases coerce,enhancers,GRanges-method
#' @exportMethod coerce
setAs(from="enhancers", to="GRanges", function(from){
  from@peaks
})
#' @rdname enhancers-class
#' @aliases show,enhancers-method
#' @param object An object of \code{"enhancers"}
#' @export
setMethod("show",
          signature(object="enhancers"),
          function(object){
            if(length(object@peaks)){
              cat("This is an object with ",
                  length(object@peaks),
                  " enhancers for ",
                  organism(object@genome))
            }else{
              cat("This is an empty object of enhancers.")
            }
          })
