snakemake@source("functions.R")
# SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)

contigs <- snakemake@input[["tcr"]] |>
    lapply(read.csv) |>
    loadContigs()

clones <- combineTCR(
    contigs,
    samples = snakemake@params[["sample_name"]],
    filterMulti = snakemake@params[["filter_chains"]]
)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat <- combineExpression(
    clones,
    seurat,
    cloneCall = "gene",
    chain = "both",
    proportion = TRUE
)

seurat@misc[["tcr_contigs"]] <- contigs
seurat@misc[["tcr_clones"]] <- clones

saveRDS(seurat, snakemake@output[["seurat"]])
