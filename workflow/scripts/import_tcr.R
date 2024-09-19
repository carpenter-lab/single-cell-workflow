snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)

# Patch for v1.12 until Bioconda has v2 with Bioconductor 3.19
checkList <- getFromNamespace("checkList", "scRepertoire")
parse10x <- getFromNamespace("pasrse10x", "scRepertoire")

loadContigs <- function(dir, format = "10X") {
    df <- checkList(dir)
    return(parse10x(df))
}

contigs <- snakemake@input[["tcr"]] |>
    lapply(read.csv) |>
    loadContigs(format = "10X")

clones <- combineTCR(
    contigs,
    samples = snakemake@params[["sample_name"]],
    filterMulti = snakemake@params[["filter_chains"]]
)

for (i in seq_along(clones)) {
    clones[[i]]$patient_id <- snakemake@params[["patient_id"]][[i]]
    clones[[i]]$condition <- snakemake@params[["condition"]][[i]]
}

assay <- list(contigs = contigs, clones = clones)

saveRDS(assay, snakemake@output[["assay"]])
