snakemake@source("functions.R")
SinkAllOutput(snakemake@log)

library(tidyverse)
library(patchwork)
library(rlang)
library(cli)

rlang::global_entrace()


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

ClusterPlot <- function(seurat, reduction, group_by = NULL, split_by = NULL, ...) {
    group_by_vals <- seurat[[group_by, drop = TRUE]]
    na_in_group_by <- any(is.na(group_by_vals) | group_by_vals == "NA")
    
    idents_keep <- if (na_in_group_by & is.null(split_by)) unique(group_by_vals[!is.na(group_by_vals)]) else NULL
    split_by_test <- split_by %||% ""
    if (group_by == split_by_test) {
        idents_keep <- NULL
        legend_pos <- "none"
    } else {
        legend_pos <- "bottom"
    }
    
    SCpubr::do_DimPlot(
        seurat,
        reduction = reduction,
        group.by = group_by,
        split.by = split_by,
        idents.keep = idents_keep,
        na.value = "grey80",
        legend.position = legend_pos,
        font.size = 10
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

DEPlot <- function(seurat, de_table, n_genes = 5, ...) {
    seurat <- tryCatch(
        SeuratObject::JoinLayers(seurat),
        error = function(e) return(seurat)
    )
    SeuratObject::LayerData(seurat, layer = "data") <- as(SeuratObject::LayerData(seurat, layer = "data"), "sparseMatrix")
    SCpubr::do_GroupwiseDEPlot(sample = seurat, de_genes = de_table, top_genes = n_genes)
}

ClusterCorrelationPlot <- function(seurat, group_by, assay, cluster, ...) {
    data <- Seurat::AverageExpression(seurat, assays = assay, group.by = group_by, layer = "scale.data")[[assay]]
    data_cor <- data |> as.matrix() |> cor()

    clus_freq <- seurat[[]] |>
            group_by(across(all_of(group_by))) |>
            summarise(n = n()) |>
            mutate(n = n / sum(n)) |>
            pull(n)

    col_annotation <- ComplexHeatmap::HeatmapAnnotation(
        "Cluster Frequency" = ComplexHeatmap::anno_barplot(
            clus_freq,
            bar_width = 0.8,
            height = unit(100, "points"),
            gp = grid::gpar(fill = "#68A691"),
            axis_param = list(gp = grid::gpar(fontsize = 20))
        ),
        annotation_name_gp = grid::gpar(fontsize = 20),
        gp = grid::gpar(fontsize = 20)
    )
    new_labels <- 0:(nrow(data_cor) - 1)

    split <- if (cluster) 4 else letters[floor(new_labels / 5) + 1]

    hm <- ComplexHeatmap::Heatmap(
        data_cor,
        col = circlize::colorRamp2(c(-1, 0, 1), c("#023C40", "white", "#B80C09")),
        na_col = "white",

        cluster_rows = cluster,
        cluster_columns = cluster,
        border = "black",
        border_gp = grid::gpar(lwd = 2),

        row_labels = new_labels,
        row_names_gp = grid::gpar(fontsize = 30),
        row_dend_side = "right",
        row_dend_width = unit(60, "points"),
        row_names_side = if_else(cluster, "right", "left"),
        row_split = split,
        row_title = " ",
        row_gap = unit(10, "points"),

        column_labels = new_labels,
        column_names_rot = 50,
        column_names_gp = grid::gpar(fontsize = 30),
        column_dend_height = unit(60, "points"),
        column_names_side = if_else(cluster, "top", "bottom"),
        column_split = split,
        column_title = " ",
        column_gap = unit(10, "points"),

        top_annotation = col_annotation,
        heatmap_legend_param = list(
            title = "Pearson's ρ",
            legend_height = unit(500, "points"),
            border = "black",
            labels_gp = grid::gpar(fontsize = 20),
            title_gp = grid::gpar(fontsize = 30),
            title_position = "leftcenter-rot",
            grid_width = unit(20, "points")
        ),

        width = unit(1000, "points"),
        height = unit(1000, "points")
    )
    return(hm)
}

BarPlot <- function(seurat, group_by = NULL, split_by = NULL, ...) {
    legend_pos <- ifelse(is.null(split_by), "none", "bottom")
    if (!is.null(split_by)) {
        new_group_by <- split_by
        split_by <- group_by
        group_by <- new_group_by
    }
    plot <- SCpubr::do_BarPlot(seurat, group.by = group_by, split.by = split_by, legend.position = legend_pos, position = "stack")
    if (!is.null(split_by)) {
        plot <- plot + SCpubr::do_BarPlot(seurat, group.by = group_by, split.by = split_by, legend.position = legend_pos, position = "fill")
    }
    return(plot)
}

PlotMethod <- function(type) {
    fun <- switch(
        type,
        qc = QCPlot,
        cluster = ClusterPlot,
        cell_type = MultiCellTypeAnnDimPlot,
        dot = DotPlot,
        heatmap = Heatmap,
        de_plot = DEPlot,
        cluster_cor = ClusterCorrelationPlot,
        bar = BarPlot
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
    valid_params[["de_table"]] <- de
}

plot <- PlotMethod(named_params[["plot"]])(seurat, !!!valid_params)

plot_title_ann <- plot_annotation(title = valid_params[["title"]], subtitle = valid_params[["subtitle"]])

if (!("ComplexHeatmap" %in% class(plot) || "ListOfComplexHeatmaps" %in% class(plot))) plot <- plot + plot_title_ann

wide_plot <- any(named_params[["plot"]] %in% c("de_plot", "heatmap"), "split_by" %in% named_params)

if (wide_plot) {
    width <- 16
    height <- 14
} else {
    width <- 14
    height <- 14
}

for (f in snakemake@output) {
    if ("HeatmapList" %in% class(plot)) {
        Cairo::Cairo(width = unit(width * 100, "points"), height = unit(height * 100, "points"), dpi = 72, file = f, type = tools::file_ext(f))
        ComplexHeatmap::draw(plot, heatmap_legend_side = if_else(named_params[["cluster"]], "left", "right"))
        graphics.off()
    } else {
        ggsave(
            filename = f,
            plot = plot,
            width = width,
            device = if (tools::file_ext(f) == "svg") Cairo::CairoSVG else NULL,
            height = height,
            dpi = if (tools::file_ext(f) == "png") 150 else 600
        )
    }
}
