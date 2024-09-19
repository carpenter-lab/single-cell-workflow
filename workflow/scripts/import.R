snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(Seurat)
library(scRepertoire)
library(glue)
library(patchwork)
library(data.table)
library(foreach)
library(rlang)
library(cli)

foreach2 <- function(..., .combine = "merge") {
    do.call(foreach::foreach, list2(..., .combine = .combine))
}
h5_assay_name <- snakemake@params[["h5_assay_name"]]
seurat_object_name <- switch(h5_assay_name, "Gene Expression" = "RNA", "Antibody Capture" = "ADT")
snakemake@params[["h5_assay_name"]] <- NULL

samples <- list_c(list(
    snakemake@input[names(snakemake@input) != ""],
    snakemake@params[names(snakemake@params) != ""]
))

sample_list_names <- c("h5", "patient_id", "condition", "sample_name")
names(sample_list_names) <- rep_along("*", sample_list_names)

if (any(!(sample_list_names %in% names(samples)))) {
    cli_abort(c(
        "x" = "Missing at least one sample attribute!",
        "Include all the following attributes:",
        sample_list_names
    ))
}

assay <- foreach2(!!!samples) %do% {
    data <- Read10X_h5(h5)
    object <- CreateAssay5Object(data[[h5_assay_name]])

    object <- CreateSeuratObject(object, assay = seurat_object_name)
    object <- RenameCells(object, add.cell.id = sample_name)
    object[["condition"]] <- condition
    object[["patient_id"]] <- patient_id
    object[["original_sample_name"]] <- sample_name
    object
}

assay <- JoinLayers(assay)
assay[[seurat_object_name]] <- split(assay[[seurat_object_name]], assay$original_sample_name)

saveRDS(assay, file = snakemake@output[["assay"]])