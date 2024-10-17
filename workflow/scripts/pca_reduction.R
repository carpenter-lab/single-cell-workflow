snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)

options(future.globals.maxSize = 4 * 1024^3)

RunPCAWithBPCells <- function(
    object,
    assay = NULL,
    features = NULL,
    npcs = 50,
    weight.by.var = TRUE,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 42,
    layer = "scale.data",
    verbose = FALSE
) {
    if (!is.null(x = seed.use)) set.seed(seed = seed.use)

    assay <- assay %||% DefaultAssay(object)
    matrix_use <- BPCells::write_matrix_dir(LayerData(object, layer = layer, assay = assay), tempfile("mat"))
    svd <- BPCells::svds(t(matrix_use), k = npcs)
    cell_embeddings <- if (weight.by.var) BPCells::multiply_cols(svd$u, svd$d) else svd$u

    feature_loadings <- svd$v
    rownames(feature_loadings) <- rownames(matrix_use)
    colnames(feature_loadings) <- paste0("PC_", 1:50)

    rownames(x = cell_embeddings) <- colnames(x = matrix_use)
    colnames(x = cell_embeddings) <- colnames(x = feature_loadings)

    sdev <- svd$d / sqrt(max(1, nrow(x = matrix_use) - 1))
    total_variance <- sum(BPCells::matrix_stats(matrix_use, row_stats = "variance")$row_stats["variance",])

    reduction <- CreateDimReducObject(
        embeddings = cell_embeddings,
        loadings = feature_loadings,
        assay = assay,
        stdev = sdev,
        key = reduction.key,
        misc = list(total.variance = total_variance)
    )
    object[[reduction.name]] <- reduction
    object <- LogSeuratCommand(object)
    return(object)
}

seurat <- readRDS(snakemake@input[["seurat"]])

var_features <- VariableFeatures(seurat)
unwanted_genes <- grep(pattern = "^TR[ABDG][VDJ]", x = var_features, value = TRUE)
new_var_features <- var_features[!(var_features %in% unwanted_genes)]
VariableFeatures(seurat) <- new_var_features

data_backend <- class(LayerData(seurat, "scale.data")) |> attr("package")

RunPCA <- ifelse(data_backend == "BPCells", RunPCAWithBPCells, Seurat::RunPCA)

seurat <- RunPCA(seurat, verbose = FALSE)

saveRDS(seurat, file = snakemake@output[["seurat"]])