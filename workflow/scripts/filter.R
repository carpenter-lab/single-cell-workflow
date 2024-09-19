snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat <- subset(seurat, subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    nCount_RNA < 18000 &
    percent.mt < 10)

saveRDS(seurat, file = snakemake@output[["seurat"]])