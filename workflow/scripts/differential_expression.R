snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(future)
library(data.table)

options(future.globals.maxSize = 4 * 1024^3)

plan(multicore, workers = snakemake@threads)

seurat <- readRDS(snakemake@input[["seurat"]])

DefaultAssay(seurat) <- snakemake@params[["assay"]]
Idents(seurat) <- snakemake@params[["idents"]]

seurat <- tryCatch(
    JoinLayers(seurat),
    error = function(e) return(seurat)
)

if (snakemake@params[["assay"]] == "SCT") {
    seurat <- PrepSCTFindMarkers(seurat, verbose = FALSE)
}

params <- c(snakemake@wildcards, snakemake@params)
named_params <- params[names(params) != ""]
valid_params <- named_params[named_params != "None"]

if ("ident2" %in% names(valid_params)) {
    if (!("ident1" %in% names(valid_params))) {
        cli::cli_abort(c("x" = "ident1 must be provided if ident2 is provided"))
    }
    markers <- FindMarkers(
        seurat,
        test.use = valid_params[["test_use"]],
        ident.1 = valid_params[["ident1"]],
        ident.2 = valid_params[["ident2"]],
        verbose = FALSE
    )
} else if ("ident1" %in% names(valid_params)) {
    markers <- FindMarkers(
        seurat,
        test.use = valid_params[["test_use"]],
        ident.1 = valid_params[["ident1"]],
        verbose = FALSE
    )
} else {
    markers <- FindAllMarkers(
        seurat,
        test.use = snakemake@params[["test_use"]],
        verbose = FALSE
    )
}

setDT(markers)

sep_use <- switch(
    tools::file_ext(snakemake@output[["markers"]]),
    tsv = "\t",
    csv = ",",
    "\t"
)

fwrite(markers, file = snakemake@output[["markers"]], sep = sep_use)
