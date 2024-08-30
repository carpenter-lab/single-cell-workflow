snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(Seurat)
library(scRepertoire)
library(Trex)
library(glue)
library(patchwork)
library(data.table)
library(foreach)
library(rlang)
library(cli)

foreach2 <- function(..., .combine = "merge") {
    do.call(foreach::foreach, list2(..., .combine = .combine))
}

samples <- list_c(list(
    snakemake@input[names(snakemake@input) != ""],
    snakemake@params[names(snakemake@params) != ""]
))

sample_list_names <- c("gex", "tcr", "patient_id", "condition", "sample_name")
names(sample_list_names) <- rep_along("*", sample_list_names)

if (any(!(sample_list_names %in% names(samples)))) {
    cli_abort(c(
        "x" = "Missing at least one sample attribute!",
        "Include all the following attributes:",
        sample_list_names
    ))
}

seurat <- foreach2(!!!samples) %do% {
    data <- Read10X_h5(gex)
    object <- CreateSeuratObject(data[["Gene Expression"]])
    object[["ADT"]] <- CreateAssay5Object(data[["Antibody Capture"]])
    tcr <- read.csv(tcr) |>
        loadContigs() |>
        combineTCR(samples = sample_name, removeNA = TRUE, filterMulti = TRUE)

    object <- RenameCells(object, add.cell.id = sample_name)

    object <- combineExpression(
        tcr,
        object,
        cloneCall = 'gene',
        chain = "both",
        proportion = TRUE
    )
    object[["condition"]] <- condition
    object[["patient_id"]] <- patient_id

    object
}

saveRDS(seurat, file = snakemake@output[["seurat"]])