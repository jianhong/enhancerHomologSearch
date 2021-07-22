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
    stopifnot("The length of subject must be less than 3" = length(subject)<3)
    null <- lapply(subject, FUN = function(.ele){
      stopifnot("subject must be an list of object of Enhancers" =
                  is(.ele, "Enhancers"))
    })
  }else{
    stopifnot("subject must be an object of Enhancers" =
                is(subject, "Enhancers"))
  }
  stopifnot("query must be an object of DNAStringSet" =
              is(query, "DNAStringSet"))
  stopifnot("The length of query must be 1" = length(query)==1)
}
