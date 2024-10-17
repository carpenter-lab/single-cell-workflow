snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)

options(future.globals.maxSize = 4 * 1024^3)

seurat <- readRDS(snakemake@input[["seurat"]])
method <- switch(
    snakemake@params[["method"]],
    "harmony" = HarmonyIntegration,
    "cca" = CCAIntegration,
    "rpca" = FastRPCAIntegration
)

seurat <- IntegrateLayers(
    seurat,
    method = method,
    new.reduction = snakemake@params[["method"]],
    dims = 1:20,
    verbose = FALSE
)

saveRDS(seurat, file = snakemake@output[["seurat"]])