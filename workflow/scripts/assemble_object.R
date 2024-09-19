snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)
library(cli)

inputs <- snakemake@input[names(snakemake@input) != ""]

assays <- lapply(inputs, readRDS)

first_assay <- NULL

if ("RNA" %in% names(assays)) {
    seurat <- assays[["RNA"]]
    first_assay <- "RNA"
} else if ("ADT" %in% names(assays)) {
    seurat <- assays[["ADT"]]
    first_assay <- "ADT"
} else {
    cli_abort("Must have RNA or ADT Assay present to use this script!")
}

if (first_assay == "RNA") {
    if ("ADT" %in% names(assays)) seurat[["ADT"]] <- assays[["ADT"]][["ADT"]]
}

if ("TCR" %in% names(assays)) {
    seurat <- combineExpression(
        assays[["TCR"]][["clones"]],
        seurat,
        cloneCall = "gene",
        chain = "both",
        proportion = TRUE
    )

    seurat@misc[["tcr_contigs"]] <- assays[["TCR"]][["contigs"]]
    seurat@misc[["tcr_clones"]] <- assays[["TCR"]][["clones"]]
}

saveRDS(seurat, snakemake@output[["seurat"]])