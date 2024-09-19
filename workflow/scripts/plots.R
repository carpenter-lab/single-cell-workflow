snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(patchwork)
library(rlang)
library(cli)
options(error = rlang::entrace)


withr::with_options(
    list(rlang_backtrace_on_error = "none"),
    if (length(snakemake@output) > 1) {
        cli_abort(c(
            "This script only accepts 1 output file.",
            "x" = "{length(snakemake@output)} output files detected",
            "i" = "Please set `output` in the Snakemake rule to only one file."
        ))
    }
)


QCPlot <- function(seurat, ...) {
    p1 <- Seurat::FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
    p2 <- Seurat::FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    p3 <- map2(
        c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        list(200, 18000, 10),
        ~SCpubr::do_ViolinPlot(seurat, features = .x, y_cut = .y)
    ) |> wrap_plots(nrow = 1)
    p3[[1]] <- p3[[1]] +
        geom_hline(
            yintercept = 4000,
            linetype = "longdash",
            colour = "black",
            linewidth = 1,
            na.rm = TRUE
        )

    (p1 + p2) / p3 & Seurat::NoLegend()
}

ClusterPlot <- function(seurat, reduction, group_by = NULL, ...) {
    Seurat::DimPlot(seurat, reduction = reduction, group.by = group_by) & theme(aspect.ratio = 1)
}

SplitDimPlot <- function(seurat, reduction = "umap", group_by = NULL, split_by = "condition", ...) {
    shuffle <- TRUE
    idents_keep <- NULL

    if (!is.null(group_by)) {
        split_levels <- unique(seurat[[group_by, drop = TRUE]])
        idents_keep <- c(split_levels[is.na(split_levels)], split_levels[!is.na(split_levels)])
        shuffle <- !any(is.na(idents_keep))
        idents_keep <- if (!all(is.na(idents_keep))) NULL
    }

    Seurat::DimPlot(
        seurat,
        reduction = reduction,
        group.by = group_by,
        split.by = split_by
    ) &
        theme(aspect.ratio = 1)
}

MultiCellTypeAnnDimPlot <- function(seurat, reduction, split_by = NULL, group_by = NULL, ...) {
    ann_level1 <- split_by
    ann_level2 <- group_by

    ann_levels <- unique(seurat[[ann_level1, drop = TRUE]])
    plots <- list()
    SeuratObject::Idents(seurat) <- ann_level1

    for (ann in ann_levels) {
        seurat_mod <- seurat
        seurat_mod[[]][SeuratObject::WhichCells(seurat_mod, idents = ann, invert = TRUE), ann_level2] <- NA

        plots[[ann]] <- Seurat::DimPlot(
            seurat_mod,
            reduction = reduction,
            group.by = ann_level2,
            shuffle = TRUE
        ) +
            theme(aspect.ratio = 1) +
            ggtitle(ann)
    }
    wrap_plots(plots) &
        Seurat::NoAxes() &
        theme(
            legend.key.size = unit(1, 'pt'),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 5)
        )
}

DotPlot <- function(seurat, features, group_by = NULL, assay = NULL, ...) {
    SCpubr::do_DotPlot(
        seurat,
        features = features,
        flip = TRUE,
        cluster = TRUE,
        group.by = group_by,
        assay = assay
    )
}

Heatmap <- function(seurat, group_by, features = NULL, ...) {
    Seurat::DoHeatmap(seurat, group.by = group_by, features = features)
}

PlotMethod <- function(type) {
    fun <- switch(
        type,
        qc = QCPlot,
        cluster = ClusterPlot,
        split = SplitDimPlot,
        cell_type = MultiCellTypeAnnDimPlot,
        dot = DotPlot,
        heatmap = Heatmap
    )
    function(...) do.call(fun, list2(...))
}

seurat <- readRDS(snakemake@input[["seurat"]])

params <- c(snakemake@wildcards, snakemake@params)
named_params <- params[names(params) != ""]
valid_params <- named_params[named_params != "None"]

if ("de" %in% names(snakemake@input)) {
    de <- read_tsv(snakemake@input[["de"]])
    features <- de |>
        group_by(cluster) |>
        slice_max(order_by = avg_log2FC, n = 30) |>
        pull(gene)
    valid_params[["features"]] <- features
}

plot <- PlotMethod(named_params[["plot"]])(seurat, !!!valid_params) +
    plot_annotation(
        title = valid_params[["title"]],
        subtitle = valid_params[["subtitle"]]
    )

ggsave(
    filename = snakemake@output[[1]],
    plot = plot,
    width = 14, height = 10, dpi = 600
)
