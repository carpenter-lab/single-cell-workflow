snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)
library(Trex)

options(future.globals.maxSize = 4 * 1024^3)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat <- quietTCRgenes(seurat)
seurat <- RunPCA(seurat, verbose = FALSE)

saveRDS(seurat, file = snakemake@output[["seurat"]])