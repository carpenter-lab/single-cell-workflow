snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)

options(future.globals.maxSize = 4 * 1024^3)

seurat <- readRDS(snakemake@input[["seurat"]])

if ("CTgene" %in% colnames(seurat[[]])) {
    seurat <- Trex::quietTCRgenes(seurat)
}

seurat <- RunPCA(seurat, verbose = FALSE)

saveRDS(seurat, file = snakemake@output[["seurat"]])