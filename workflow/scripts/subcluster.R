snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(Seurat)
library(rlang)


SubsetByTCR <- function(seurat, ...) {
    seurat[["t_cell_receptor"]] <- !is.na(seurat$CTaa)
    subset(seurat, cells = WhichCells(seurat, expression = t_cell_receptor))
}

SubsetByIdent <- function(seurat, ident_use, idents_keep, ...) {
    print("here")
    Idents(seurat) <- ident_use
    subset(seurat, cells = WhichCells(seurat, idents = idents_keep))
}

seurat <- readRDS(snakemake@input[["seurat"]])

params <- c(snakemake@wildcards, snakemake@params)
named_params <- params[names(params) != ""]
valid_params <- named_params[named_params != "None"]

SubsetMethod <- function(seurat, ...) {
    dots <- list2(...)
    if (!is_missing("ident_use") && !is.null(dots[["ident_use"]])) {
        check_required("idents_keep")
        f <- SubsetByIdent
    } else {
        f <- SubsetByTCR
    }
    dots[["seurat"]] <- seurat
    do.call(f, dots)
}

seurat <- SubsetMethod(seurat, !!!valid_params)

saveRDS(seurat, file = snakemake@output[["seurat"]])