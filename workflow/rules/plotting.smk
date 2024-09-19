from common import (
    report_plot_labels,
    get_proper_clustering_output,
    get_dot_plot_features,
    PlotTitle,
    get_plot_type,
    get_category_name,
    make_plot_subtitle,
)


def qc_plot_title(wildcards):
    if wildcards["type"] == "pre":
        return "QC Plots Before Filtering"
    else:
        return "QC Plots After Filtering"


def qc_plot_input(wildcards):
    if wildcards["type"] == "pre":
        return "results/qc/prep.rds"
    else:
        return "results/seurat_objects/all_data.rds"


rule qc_plot:
    input:
        seurat=qc_plot_input,
    params:
        title=qc_plot_title,
        plot="qc",
    output:
        report(
            "results/qc/{type}_qc_plot.png",
            category="Preprocessing",
            subcategory="Quality Control",
            labels=report_plot_labels,
        ),
    log:
        "logs/qc/{type}_plot.log",
    conda:
        "../envs/plotting.yml"
    script:
        "../scripts/plots.R"


use rule qc_plot as umap_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
    params:
        title=PlotTitle(get_plot_type).make_title,
        subtitle=make_plot_subtitle,
        plot=get_plot_type,
        reduction="{assay}_{reduction_use}_umap",
    output:
        report(
            "results/clustering/{subset}/plots/{assay}/{reduction_use}/{group_by}_split_{split_by}.png",
            category=get_category_name,
            subcategory="Clustering",
            labels=report_plot_labels,
        ),
    log:
        "logs/umap_plots/{subset}/{assay}/{reduction_use}/{group_by}_split_{split_by}.log",


use rule qc_plot as dot_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
    params:
        title=PlotTitle("dot").make_title,
        subtitle=make_plot_subtitle,
        plot="dot",
        features=get_dot_plot_features(config),
    output:
        report(
            "results/de/{subset}/plots/{assay}/dot_plot_by_{group_by}.png",
            category=get_category_name,
            subcategory="Differential Expression",
            labels=report_plot_labels,
        ),
    log:
        "logs/de/{subset}/plots/{assay}/dot_plot_by_{group_by}.log",


use rule qc_plot as heatmap with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        de="results/de/{subset}/{assay}/{group_by}_markers.tsv",
    params:
        title=PlotTitle("heatmap").make_title,
        subtitle=make_plot_subtitle,
        plot="heatmap",
    output:
        report(
            "results/de/{subset}/plots/{assay}/heatmap_by_{group_by}.png",
            category=get_category_name,
            subcategory="Differential Expression",
            labels=report_plot_labels,
        ),
    resources:
        mem="5GB",
    log:
        "logs/de/{subset}/plots/{assay}/heatmap_by_{group_by}.log",
