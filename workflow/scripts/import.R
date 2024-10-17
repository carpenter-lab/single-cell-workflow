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
library(BPCells)

foreach2 <- function(..., .combine = "merge") {
    do.call(foreach::foreach, list2(..., .combine = .combine))
}

h5_assay_name <- snakemake@params[["h5_assay_name"]]
seurat_object_name <- switch(h5_assay_name, "Gene Expression" = "RNA", "Antibody Capture" = "ADT")
snakemake@params[["h5_assay_name"]] <- NULL

if ("use_bpcells" %in% names(snakemake@params) && snakemake@params[["use_bpcells"]]) {
    if (length(snakemake@output[["matrix_dir"]]) != 1) {
        cli_abort("If use_bpcells is true, a single directory must be provided as an input.")
    }

    if (h5_assay_name != "Gene Expression") {
        cli_abort("Using BPCells for assays other than RNA or ATAC is not reccomended.")
    }
} else {
    # Set to empty string to prevent errors later. 
    # since use_bpcells is FALSE, the directory won't be used
    snakemake@output[["matrix_dir"]] <- ""
}
use_bpcells <- snakemake@params[["use_bpcells"]] %||% FALSE
snakemake@params[["use_bpcells"]] <- NULL

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

Read10X <- function(h5, use_bpcells = FALSE, matrix_dir = NULL, assay = NULL) {
    data <- Read10X_h5(h5)[[assay]]
    if (use_bpcells) {
        BPCells::write_matrix_dir(mat = data, dir = matrix_dir)
        data <- BPCells::open_matrix_dir(dir = matrix_dir)
    }
    return(data)
}

assay <- foreach2(!!!samples) %do% {
    matrix_dir_use <- file.path(snakemake@output[["matrix_dir"]], sample_name)
    data <- Read10X(h5 = h5, use_bpcells = use_bpcells, matrix_dir = matrix_dir_use, assay = h5_assay_name)
    object <- CreateAssay5Object(data)

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