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
splitStrByBlock <- function(s, block, FUN, FUN.VALUE=character(1L), ...){
  s <- split(s, rep(seq.int(ceiling(length(s)/block)),
                    each=block)[seq_along(s)])
  vapply(s, FUN=FUN, FUN.VALUE = FUN.VALUE, ...)
}
write.phylip_html <- function(x, filepath, block=80){
  x <- as.character(x)
  x100 <- strsplit(x, split="")
  n <- length(x100[[1]])
  N <- ceiling(n/block)
  x100 <- lapply(x100, function(.ele) {
    splitStrByBlock(.ele, block=block, FUN=paste, collapse="")
  })
  blank <- c("<br/>", "<br/>")
  x100 <- do.call(rbind, x100)
  x100 <- apply(x100, 2, FUN=function(.ele){
    paste("<p>", c(.ele, blank), "</p>")
  })
  x100_n <- paste("<p>", c(names(x), blank), "</p>")
  x100_n <- rep(c("<div>", x100_n, "</div>"), N)
  counter <- lapply(unique(c(seq.int(N-1)*block, n)), function(.ele){
    paste("<p>",
          c(.ele, rep("<br/>", length(x)-1), blank),
          "</p>")
  })
  counter <- unlist(counter)
  html <- c("<!DOCTYPE html><html><head><title>",
            basename(filepath),
            "</title>",
            "<style>",
            "#align { margin-left: calc(50% - 50vw); width: 100vw; font-family: monospace;}",
            "#seqnames,#sequences,#counts { float: left}",
            "#seqnames { max-width: 20%; }",
            "#sequences {width: auto; min-width: 700px; margin-left: 20px;}",
            "#counts { width: auto; align: 'left'}",
            ".space { margin-right: 12px;}",
            "p {margin: 0px;}",
            "</style>",
            '<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>',
            "<script>",
            '$(document).ready(function(){',
            "var currentMousePos = { x: -1, y: -1};",
            '$("#sequences").html(function(i,el) {',
            "return el.replace(/((A|C|G|T|N|-){10})/g, '$1<span class=\"space\"></span>')",
            ".replace(/A/g, '<span style=\"color:green\">A</span>')",
            ".replace(/C/g, '<span style=\"color:orange\">C</span>')",
            ".replace(/G/g, '<span style=\"color:red\">G</span>')",
            ".replace(/T/g, '<span style=\"color:blue\">T</span>')",
            ".replace(/U/g, '<span style=\"color:blue\">U</span>');",
            '});',
            "});",
            "</script>",
            "</head><body>",
            "<h1>Enhancer homologs multiple alignments</h1>",
            "<div id='align'>",
            "<div id='seqnames'>", x100_n, "</div>",
            "<div id='sequences'>", x100, "</div>",
            "<div id='counts'>", counter, "</div>",
            "</div></body></html>")
  writeLines(html, filepath)
}
#' @importFrom Biostrings consensusMatrix DNAStringSet write.phylip
#' matchPattern alphabetFrequency injectHardMask
#' @importFrom IRanges IRangesList Views
#' @importFrom utils write.csv
saveAln <- function (x, filepath, motifConsensus, format) {
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
    involved_motifs <- sort(mt)
    mt <- reduce(mt)
    if(length(mt)){
      m <- Views(cons_m, start = start(mt), end = end(mt))
      m <- injectHardMask(m, letter="-")
      m <- DNAStringSet(m)
      names(m) <- "motifConsensus"
      x@unmasked <- c(x@unmasked, m)
      write.csv(involved_motifs, file = paste0(sub(paste0(".",format,"$"),
                                                   "", filepath),
                                               ".motifConsensus.info.csv"),
                row.names=FALSE)
    }
  }
  if(format=="txt"){
    write.phylip(x, filepath = filepath)
  }else{
    write.phylip_html(x, filepath = filepath)
  }
}
#' output alignments
#' @description Save enhancer homologs to file in phylip format.
#' @param al output of \link{alignment}.
#' @param output_folder output folder.
#' @param motifConsensus Transcription factor binding consensus.
#' @param format The format of output files.
#'  Available formats are 'txt' and 'html'. Default is 'txt'.
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
saveAlignments <- function(al,
                           output_folder=tempdir(),
                           motifConsensus=NULL,
                           format=c("txt", "html")){
  checkMotifCons(motifConsensus)
  format <- match.arg(format)
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  if(length(names(al))!=length(al)){
    names(al) <- paste0("aln", seq_along(al), ".phylip.", format)
  }else{
    names(al) <- paste0(names(al), ".phylip.", format)
  }
  null <- mapply(saveAln, x=al,
                 filepath = file.path(output_folder, names(al)),
                 MoreArgs=list(motifConsensus = motifConsensus,
                               format=format))
  return(invisible(null))
}
