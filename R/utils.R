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
