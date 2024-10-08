snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(future)
library(rlang)

options(
    future.globals.maxSize = 4 * 1024^3,
    future.rng.onMisuse = "ignore", # silences an erronious RNG Parrelezation error
    Seurat.warn.umap.uwot = FALSE # silences the warning that Seurat uses UWOT instead of Python
)

RunDimensionReductionAndCluster <- function(object, assay = NULL, reduction = NULL) {
    current_assay <- SeuratObject::DefaultAssay(object)
    on.exit(expr = SeuratObject::DefaultAssay(object) <- current_assay)

    graph_types <- c("nn", "snn")

    SeuratObject::DefaultAssay(object) <- assay <- assay %||% SeuratObject::DefaultAssay(object)
    object <- Seurat::RunUMAP(
        object,
        reduction = reduction,
        dims = 1:30,
        verbose = FALSE,
        umap.method = "uwot",
        reduction.name = glue::glue("{assay}_{reduction}_umap"),
    )

    object <- Seurat::FindNeighbors(
        object,
        reduction = reduction,
        dims = 1:30,
        verbose = FALSE,
        graph.name = glue::glue("{assay}_{reduction}_{graph_types}")
    )
    object <- Seurat::FindClusters(
        object,
        verbose = FALSE,
        graph.name = glue::glue("{assay}_{reduction}_snn")
    )
    object[[glue::glue("{assay}_{reduction}_clusters")]] <- object[["seurat_clusters"]]
    object[["seurat_clusters"]] <- NULL
    return(object)
}

plan(multicore, workers = snakemake@threads)

seurat <- readRDS(snakemake@input[["seurat"]])

for (reduction in snakemake@params[["reductions"]]) {
    seurat <- RunDimensionReductionAndCluster(
        seurat,
        assay = snakemake@params[["assay"]], reduction = reduction
    )
}

saveRDS(seurat, file = snakemake@output[["seurat"]])