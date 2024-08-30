snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(future)
library(Azimuth)
library(glue)
library(stringr)

options(future.globals.maxSize = 16 * 1024^3)

plan(multicore, workers = snakemake@threads)

seurat <- readRDS(snakemake@input[["seurat"]])
reference_use <- snakemake@params[["reference"]]

azimuth_references <- SeuratData::AvailableData()$Dataset[grepl("*ref", SeuratData::AvailableData()$Dataset)]
if (!(reference_use %in% azimuth_references)) {
    names(azimuth_references) <-  rep_along("*", azimuth_references)
    cli::cli_abort(c(
        "x" = "{reference_use} is not a valid Aziumuth reference!",
        " " = "Possible references include:",
        azimuth_references
    ))
}

if (!require(glue("{reference_use}.SeuratData"), character.only = TRUE)) {
    install.packages(glue("{reference_use}.SeuratData"), repos = "https://seurat.nygenome.org/")
}

seurat <- RunAzimuth(
    seurat,
    reference = reference_use,
    umap.name = glue("azimuth_{reference_use}_umap") |> str_replace("\\.", "_"),
    verbose = FALSE,
    k.weight = 100,
    n.trees = 40,
    mapping.score.k = 100
)

saveRDS(seurat, file = snakemake@output[["seurat"]])