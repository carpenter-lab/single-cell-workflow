snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(scRepertoire)

options(future.globals.maxSize = 4 * 1024^3)

seurat <- readRDS(snakemake@input[["seurat"]])

if (snakemake@params[["method"]] == "sctransform") {
    seurat <- SCTransform(seurat, verbose = FALSE)
} else {
    seurat <- NormalizeData(marrow)
    seurat <- FindVariableFeatures(seurat)
}

seurat <- CellCycleScoring(
    seurat,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    s.features = cc.genes.updated.2019$s.genes,
    search = TRUE,
    nbin = 12
)

seurat[["cell_cycle_diff"]] <- seurat[["S.Score"]] - seurat[["G2M.Score"]]

if (length(snakemake@params[["vars_to_regress"]] > 0)) {
    if (snakemake@params[["method"]] == "sctransform") {
        seurat <- SCTransform(
            seurat,
            verbose = FALSE,
            vars.to.regress = snakemake@params[["vars_to_regress"]]
        )
    } else {
        seurat <- ScaleData(seurat, snakemake@params[["vars_to_regress"]])
    }
}

if ("ADT" %in% Assays(seurat)) {
     seurat <- NormalizeData(
         seurat,
         assay = "ADT",
         normalization.method = "CLR",
         verbose = FALSE
     )
    seurat <- ScaleData(seurat, assay = "ADT", verbose = FALSE)
    seurat <- FindVariableFeatures(seurat, assay = "ADT", verbose = FALSE)
}

saveRDS(seurat, file = snakemake@output[["seurat"]])
