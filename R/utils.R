swapList <- function (x){
  stopifnot(is.list(x))
  null <- sapply(x, function(.ele) {
    stopifnot(is.list(.ele) || is(.ele, "CompressedIRangesList"))
  })
  levelsA <- names(x)
  levelsB <- unique(unlist(sapply(x, names, simplify = FALSE)))
  if(length(levelsB)==0){
    stopifnot(all(lengths(x)==length(x[[1]])))
    levelsB <- seq.int(length(x[[1]]))
  }
  y <- as.list(levelsB)
  names(y) <- levelsB
  for (.lB in levelsB) {
    y[[.lB]] <- list()
  }
  for (.lA in levelsA) {
    for (.lB in levelsB) {
      y[[.lB]][[.lA]] <- x[[.lA]][[.lB]]
    }
  }
  y
}


#' @importFrom GenomeInfoDb mapGenomeBuilds
guessAssembly <- function(genome){
  # hard coding, maybe an issue.
  assembly <- genome@metadata$genome
  assembly <- mapGenomeBuilds(assembly, style="Ensembl")
  assembly_ensemblID <- assembly[nrow(assembly), "ensemblID"]
  assembly_ensemblID <- sub("\\.p\\d+$", "", assembly_ensemblID)
  assembly_ucscID <- assembly[nrow(assembly), "ucscID"]
  assembly <- unique(c(assembly_ensemblID, assembly_ucscID))
  return(assembly)
}

# help function to check the query and subject parameters for alignment
# functions
checkQuerySubject <- function(query, subject, subjectIsList=FALSE){
  if(subjectIsList){
    stopifnot("subject must be an list of object of Enhancers" =
                is(subject, "list"))
    stopifnot("The length of subject must be more than 0" = length(subject)>0)
    #stopifnot("The length of subject must be less than 3" = length(subject)<3)
    null <- lapply(subject, FUN = function(.ele){
      stopifnot("subject must be an list of object of Enhancers" =
                  is(.ele, "Enhancers"))
      stopifnot("Enhancers are empty" = length(peaks(.ele))>0)
    })
  }else{
    stopifnot("subject must be an object of Enhancers" =
                is(subject, "Enhancers"))
    stopifnot("Enhancers are empty" = length(peaks(subject))>0)
  }
  stopifnot("query must be an object of DNAStringSet" =
              is(query, "DNAStringSet"))
  stopifnot("The length of query must be 1" = length(query)==1)
}

# help function to check substitutionMatrix
checkSubstitutionMatrix <- function(substitutionMatrix=c("iub", "clustalw"),
                                    method=c("Muscle", "ClustalW")){
  method <- match.arg(method)
  if(any(is.na(substitutionMatrix))){
    stop("substitutionMatrix can not be NA")
  }
  if(is.null(substitutionMatrix)){
    return(substitutionMatrix)
  }
  if(method=="Muscle" &&is.character(substitutionMatrix)){
    return(NULL)
  }
  if(is.character(substitutionMatrix)){
    return(match.arg(substitutionMatrix))
  }
  if(is.matrix(substitutionMatrix)){
    seqn <- c('A', 'C', 'G', 'T')
    stopifnot("rownames of substitutionMatrix must be c('A', 'C', 'G', 'T')" =
                all(seqn %in% rownames(substitutionMatrix)))
    stopifnot("colnames of substitutionMatrix must be c('A', 'C', 'G', 'T')" =
                all(seqn %in% colnames(substitutionMatrix)))
    if(method=="Muscle"){
      headerNames <- c("A", "C", "D", "E", "F",
                       "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R",
                       "S", "T", "V", "W", "Y")
    }else{
      headerNames <- c("A", "R", "N", "D", "C",
                       "Q", "E", "G", "H", "I",
                       "L", "K", "M", "F", "P",
                       "S", "T", "W", "Y", "V",
                       "B", "Z", "X", "*")
    }
    auxMat <- matrix(0, length(headerNames), length(headerNames))
    rownames(auxMat) <- headerNames
    colnames(auxMat) <- headerNames
    auxMat[seqn, seqn] <- substitutionMatrix
    substitutionMatrix <- auxMat
  }
  return(substitutionMatrix)
}

# help function to check PWMs
checkPWMs <- function(PWMs){
  stopifnot("PWMs must be an object of PWMatrixList or PFMatrixList"=
              inherits(PWMs, c("PWMatrixList", "PFMatrixList")))
}
# help function to check TFBP bar
checkTFBPs <- function(x){
  stopifnot("Input of TFBPscore format is not correct"=is.logical(x))
  stopifnot("Input of TFBPscore format is not correct"=
              length(names(x))==length(x))
}

# help function to prepare parameters for ClustalW
prepareClustalWparams <- function(order=c("aligned", "input"),
                                  substitutionMatrix, ...){
  params <- list(...)
  params[["inputSeqIsFileFlag"]] <- FALSE # character from DNAStringSet
  params[["outorder"]] <- match.arg(order)
  # possible value:c("gcg", "gde", "pir", "phylip", "nexus", "fasta", "clustal")
  # but fixed to clustal only.
  params[["output"]] <- "clustal"
  params[["case"]] <- "upper" # "lower or upper, fixed to upper
  params[["seqnosFlag"]] <- is.null(params[["seqnos"]]) # set params[["seqnos"]] to NULL

  params[["substitutionMatrixIsDefaultFlag"]] <- FALSE
  params[["substitutionMatrixIsStringFlag"]] <- FALSE
  params[["substitutionMatrix"]] <- substitutionMatrix
  if(is.null(substitutionMatrix)||is.character(substitutionMatrix)){
    params[["substitutionMatrixIsDefaultFlag"]] <- TRUE
    params[["pwdnamatrix"]] <- substitutionMatrix
    params[["substitutionMatrix"]] <- "default"
  }else{
    if(is.matrix(substitutionMatrix)){
      params[["dnamatrix"]] <- NULL
    }
  }
  params[["pwmatrix"]] <- NULL # dna alignment only
  for(k in c("options", "check", "fullhelp", "align", "pim",
             "convert", "quicktree", "negative", "endgaps",
             "nopgap", "nohgap", "novgap", "noweights",
             "profile", "sequences", "nosecstr1", "nosecstr2",
             "kimura", "tossgaps")){
    if(is.null(params[[k]])||is.na(params[[k]])){
      params[[k]] <- FALSE
    }else{
      if(!is.logical(params[[k]])){
        stop("The parameter ", k, " must be logical, \n")
      }
    }
  }

  for(k in c("ktuple", "topdiags", "window", "pairgap", "gapdist",
             "maxdiv", "helixgap", "strandgap", "loopgap",
             "terminalgap", "helixendin", "helixendout",
             "strandendin", "strandendout", "seed")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be an integer(1).")
      }
      if(length(params[[k]])!=1 && params[[k]]!=round(params[[k]])){
        stop("The parameter ", k, " should be an integer(1).")
      }
    }
  }

  for(k in c("pwgapopen", "pwgapext", "transweight")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be numeric(1).")
      }
      if(length(params[[k]])!=1 &&!is.numeric(params[[k]])){
        stop("The parameter ", k, " should be numeric(1).")
      }
    }
  }

  for(k in c("stats", "usetree", "profile1", "profile2",
             "usetree1", "usetree2")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be a filename.")
      }
      if(length(params[[k]])!=1 &&!is.character(params[[k]])){
        stop("The parameter ", k, " should be a filename(1).")
      }
      if (!file.exists(params[[k]])){
        stop("The file for parameter ", k ," does not exist!")
      }
      if (file.info(params[[k]])$size == 0){
        stop("The file for parameter ", k ," is empty!")
      }
    }
  }

  if (!is.null(params[["hgapresidues"]])){
    ##check if hgapresidues is a string
    if (!is.character(params[["hgapresidues"]])){
      stop("The parameter hgapresidues should be a string!")
    }
  }

  options <- list("score"=c("percent", "absolute"),
                  "iteration"=c("tree", "alignment", "none"),
                  "secstrout"=c("structure", "mask", "both", "none"),
                  "outputtree"=c("nj", "phylip", "dist", "nexus"),
                  "bootlabels"=c("node", "branch"),
                  "seqnos"=c("on", "off"),
                  "seqno_range"=c("off", "on"),
                  "pwdnamatrix"=c("iub", "clustalw"))
  for(k in names(options)){
    if(!is.null(params[[k]])){
      params[[k]] <- options[[k]][match(params[[k]],
                                        options[[k]],
                                        nomatch = 1)[1]]
      if(is.na(params[[k]])){
        stop("The parameter", k, "is incorrect.\nPossible choices:",
             paste(options[[k]], collapse=", "))
      }
    }
  }
  if(!is.null(params[["range"]])){
    stopifnot(length(params[["range"]]==2))
    stopifnot("The parameter range should consist of 2 positive integers."=
                all(params[["range"]]==round(params[["range"]])))
    stopifnot("The parameter range should consist of 2 positive integers."=
                any(is.na(params[["range"]])))
    stopifnot("The parameter range should consist of 2 positive integers."=
                any(params[["range"]]<0))
    params[["range"]] <- lapply(params[["range"]], as.integer)
  }

  return(params)
}

# help function to prepare parameters for Muscle
prepareMuscleparams <- function(...){
  params <- list(...)
  params[["inputSeqIsFileFlag"]] <- FALSE # character from DNAStringSet
  for(k in c("spn", "core", "anchors")){
    if(is.null(params[[k]])||is.na(params[[k]])){
      params[[k]] <- TRUE
    }else{
      if(!is.logical(params[[k]])){
        stop("The parameter ", k, " must be logical, \n")
      }
    }
  }

  for(k in c("le", "sp", "sv", "brenner", "diags", "diags1", "diags2", "dimer",
             "noanchors", "nocore", "profile", "refine", "refinew", "spscore")){
    if(is.null(params[[k]])||is.na(params[[k]])){
      params[[k]] <- FALSE
    }else{
      if(!is.logical(params[[k]])){
        stop("The parameter ", k, " must be logical, \n")
      }
    }
  }
  if (!is.null(params[["anchors"]])){
    ##both positive
    if (params[["anchors"]] && params[["noanchors"]]){
      stop("The parameters anchors and noanchors \n",
           "can't be positive at the same time!")
    }
    ##both negative
    if (!params[["anchors"]] && !params[["noanchors"]]){
      stop("The parameters anchors and noanchors \n",
           "can't be negative at the same time!")
    }
  }
  if (!is.null(params[["core"]])){
    ##both positive
    if (params[["core"]] && params[["nocore"]]){
      stop("The parameters core and nocore \n",
           "can't be positive at the same time!")
    }
    ##both negative
    if (!params[["core"]] && !params[["nocore"]]){
      stop("The parameters core and nocore \n",
           "can't be negative at the same time!")
    }
  }

  for(k in c("anchorspacing", "diagbreak", "diaglength", "diagmargin", "hydro",
             "maxtrees", "refinewindow", "smoothwindow")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be an integer(1).")
      }
      if(length(params[[k]])!=1 && params[[k]]!=round(params[[k]])){
        stop("The parameter ", k, " should be an integer(1).")
      }
    }
  }

  for(k in c("center", "hydrofactor", "maxhours", "minbestcolscore",
             "minsmoothscore", "smoothscoreceil", "SUEFF")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be numeric(1).")
      }
      if(length(params[[k]])!=1 &&!is.numeric(params[[k]])){
        stop("The parameter ", k, " should be numeric(1).")
      }
    }
  }

  for(k in c("in1", "in2")){
    if(!is.null(params[[k]])){
      if(is.na(params[[k]])){
        stop("The parameter ", k, " should be a filename.")
      }
      if(length(params[[k]])!=1 &&!is.character(params[[k]])){
        stop("The parameter ", k, " should be a filename(1).")
      }
      if (!file.exists(params[[k]])){
        stop("The file for parameter ", k ," does not exist!")
      }
      if (file.info(params[[k]])$size == 0){
        stop("The file for parameter ", k ," is empty!")
      }
    }
  }
  options <- list("cluser1"=
                    c("upgma", "upgmamax", "upgmamin",
                      "upgmb", "neighborjoining"),
                  "cluster2"=
                    c("upgma", "upgmamax", "upgmamin",
                      "upgmb", "neighborjoining"),
                  "distance1"=
                    c("kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3", "kmer4_6"),
                  "distance2"=c("pctidkimura", "pctidlog"),
                  "objscore"=c("dp", "ps", "sp", "spf", "spm", "xp"),
                  "root1"=c("pseudo", "midlongestspan", "minavgleafdist"),
                  "root2"=c("pseudo", "midlongestspan", "minavgleafdist"),
                  "weight1"=c("none", "henikoff", "henikoffpb",
                              "gsc", "clustalw", "threeway"),
                  "weight2"=c("none", "henikoff", "henikoffpb",
                              "gsc", "clustalw", "threeway"))
  for(k in names(options)){
    if(!is.null(params[[k]])){
      params[[k]] <- match(params[[k]], options[[k]], nomatch = 1)[1]
    }
  }

  return(params)
}

#' @importFrom IRanges IRanges
#' @importFrom S4Vectors extractROWS
#' @importFrom utils head tail
parseAln <- function(aln){
  version <- aln[1]
  if (length(aln) < 3 ||
      !identical(sub("^\\s+$", "", aln[2:3]), c("", "")))
    stop("There is an invalid aln file!")

  aln <- tail(aln, -3)
  aln <- sub("^(\\S+\\s+\\S+)\\s*\\d*$", "\\1", aln)

  markupPattern <- "^(\\s|\\*|:|\\.)*$"

  markupLines <- grep(markupPattern, aln, perl=TRUE)
  alnLines <- gaps(as(markupLines, "IRanges"), start=1, end=length(aln))
  nseq <- unique(width(alnLines))

  if (length(nseq) != 1)
    stop("There are missing alignment rows!")

  aln <- extractROWS(aln, alnLines)
  spaces <- regexpr("\\s+", aln)
  ids <- substr(aln, 1L, spaces - 1L)
  nsplits <- length(aln) %/% nseq

  if (!identical(ids, rep.int(head(ids, nseq), nsplits)))
    stop("The alignment rows are out of order!")

  alns <- substr(aln, spaces + attr(spaces, "match.length"), nchar(aln))

  chrs <- structure(do.call(paste,
                            c(split(alns, rep(seq_len(nsplits),
                                              each=nseq)), sep="")),
                    names=head(ids, nseq))

  out <- new("DNAMultipleAlignment",
             unmasked=do.call("DNAStringSet", list(chrs)),
             rowmask=as(IRanges(), "NormalIRanges"),
             colmask=as(IRanges(), "NormalIRanges"))
}

# help function of msa
msa <- function(inputSeqs,
                method=c("Muscle", "ClustalW"),
                cluster=c("nj", "upgma", "upgmamax", "upgmamin", "upgmb"),
                substitutionMatrix=c("iub", "clustalw"),
                gapOpening=ifelse(method[1]=="ClustalW", 15.0, 400),
                gapExtension=ifelse(method[1]=="ClustalW", 6.66, 0),
                maxiters=ifelse(method[1]=="ClustalW", 3, 16),
                order=c("aligned", "input"),
                ...){
  method <- match.arg(method)
  cluster <- match.arg(cluster)
  if(is.character(substitutionMatrix) && !is.na(substitutionMatrix)){
    substitutionMatrix <- match.arg(substitutionMatrix)
  }
  substitutionMatrix <- checkSubstitutionMatrix(substitutionMatrix, method)
  order <- match.arg(order)
  if(method=="Muscle"){
    if(cluster=="nj") cluster <- "neighborjoining"
  }else{
    if(!cluster %in% c("nj", "upgma")){
      warning("When method is ClustalW, ",
              "the possible cluster is 'nj' or 'upgma'.",
              "cluster set to 'nj'")
      cluster <- 'nj'
    }
  }
  stopifnot(is.numeric(gapOpening) &&
              length(gapOpening)==1 &&
              !is.na(gapOpening))
  stopifnot(is.numeric(gapExtension) &&
              length(gapExtension)==1 &&
              !is.na(gapExtension))
  stopifnot(is.numeric(maxiters) &&
              length(maxiters)==1 &&
              !is.na(maxiters)[1] &&
              maxiters==round(maxiters))
  inputSeqNames <- names(inputSeqs)
  inputSeqs <- as.character(inputSeqs)
  names(inputSeqs) <- paste0("S", seq_along(inputSeqs))
  if(method=="Muscle"){
    params <- prepareMuscleparams(...)
    result <- RMuscle(inputSeqs, cluster, abs(gapOpening),
                      abs(gapExtension), maxiters, substitutionMatrix,
                      "dna", FALSE, params)
  }else{# ClustalW
    params <- prepareClustalWparams(order, substitutionMatrix, ...)
    substitutionMatrix <- params[["substitutionMatrix"]]
    result <- RClustalW(inputSeqs, cluster, abs(gapOpening),
                        abs(gapExtension), maxiters, substitutionMatrix,
                        "dna", FALSE, params)
  }
  out <- parseAln(result$msa)
  if (length(inputSeqNames) > 0){
    if (order == "aligned"){
      perm <- match(names(out@unmasked), names(inputSeqs))
      names(out@unmasked) <- inputSeqNames[perm]
    }else{
      perm <- match(names(inputSeqs), names(out@unmasked))
      out@unmasked <- out@unmasked[perm]
      names(out@unmasked) <- inputSeqNames
    }
  }else
    names(out@unmasked) <- NULL
  out
}
