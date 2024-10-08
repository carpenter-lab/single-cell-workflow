from common import (
    report_plot_labels,
    get_proper_clustering_output,
    get_dot_plot_features,
    PlotTitle,
    get_plot_type,
    get_category_name,
    make_plot_subtitle,
    valid_dict_key,
    WorkflowResults,
)

PLOT_FILE_TYPES = ["pdf", "svg"]


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
        matrix_dir="results/import/bpcells_backing" if config["preprocessing"]["use_bpcells"] else None
    params:
        title=qc_plot_title,
        plot="qc",
    output:
        report(
            "results/qc/violin/{type}_qc_plot.png",
            category="Preprocessing",
            subcategory="Quality Control",
            labels=report_plot_labels,
        ),
        expand("results/qc/violin/{{type}}_qc_plot.{ext}", ext=PLOT_FILE_TYPES),
    log:
        "logs/qc/violin/{type}_plot.log",
    conda:
        "../envs/plotting.yml"
    script:
        "../scripts/plots.R"


use rule qc_plot as umap_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        matrix_dir="results/import/bpcells_backing" if config["preprocessing"]["use_bpcells"] else None
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
        expand("results/clustering/{{subset}}/plots/{{assay}}/{{reduction_use}}/{{group_by}}_split_{{split_by}}.{ext}", ext=PLOT_FILE_TYPES)
    log:
        "logs/umap_plots/{subset}/{assay}/{reduction_use}/{group_by}_split_{split_by}.log",


use rule qc_plot as dot_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        matrix_dir="results/import/bpcells_backing" if config["preprocessing"]["use_bpcells"] else None
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
        expand("results/de/{{subset}}/plots/{{assay}}/dot_plot_by_{{group_by}}.{ext}", ext=PLOT_FILE_TYPES)
    log:
        "logs/de/{subset}/plots/{assay}/dot_plot_by_{group_by}.log",


use rule qc_plot as heatmap with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        de="results/de/{subset}/{assay}/{group_by}_markers.tsv",
        matrix_dir="results/import/bpcells_backing" if config["preprocessing"]["use_bpcells"] else None
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
        expand("results/de/{{subset}}/plots/{{assay}}/heatmap_by_{{group_by}}.{ext}", ext=PLOT_FILE_TYPES)
    resources:
        mem="15GB",
    log:
        "logs/de/{subset}/plots/{assay}/heatmap_by_{group_by}.log",


use rule qc_plot as bar_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        title=PlotTitle("bar").make_title,
        subtitle=make_plot_subtitle,
        plot="bar",
    output:
        report(
            "results/clustering/{subset}/plots/{assay}/{reduction}/barplot/{group_by}_by_{split_by}.png",
            category=get_category_name,
            subcategory="Quality Control",
            labels=report_plot_labels,
        ),
        expand(
            "results/clustering/{{subset}}/plots/{{assay}}/{{reduction}}/barplot/{{group_by}}_by_{{split_by}}.{ext}",
            ext=PLOT_FILE_TYPES,
        ),
    resources:
        mem="15GB",
    log:
        "logs/clustering/{subset}/plots/{assay}/{reduction}/barplot_{group_by}_by_{split_by}.log",


use rule qc_plot as de_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        de="results/de/{subset}/{assay}/{group_by}_markers.tsv",
        matrix_dir="results/import/bpcells_backing" if config["preprocessing"]["use_bpcells"] else None
    params:
        title=PlotTitle("de_plot").make_title,
        subtitle=make_plot_subtitle,
        plot="de_plot",
        n_genes=7
    output:
        report(
            "results/de/{subset}/plots/{assay}/de_plot_by_{group_by}.png",
            category=get_category_name,
            subcategory="Differential Expression",
            labels=report_plot_labels,
        ),
        expand("results/de/{{subset}}/plots/{{assay}}/de_plot_{{group_by}}.{ext}", ext=PLOT_FILE_TYPES)
    resources:
        mem="15GB",
    log:
        "logs/de/{subset}/plots/{assay}/de_plot_by_{group_by}.log",


use rule qc_plot as cluster_cor_plot with:
    input:
        seurat=get_proper_clustering_output(config, rules),
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        title=PlotTitle("cluster_cor").make_title,
        plot="cluster_cor",
        group_by="{assay}_{reduction_use}_clusters",
    output:
        report(
            "results/clustering/{subset}/plots/{assay}/{reduction_use}/correlation.png",
            category=get_category_name,
            subcategory="Clustering",
            labels=report_plot_labels,
        ),
        expand(
            "results/clustering/{{subset}}/plots/{{assay}}/{{reduction_use}}/correlation.{ext}",
            ext=PLOT_FILE_TYPES,
        ),
    resources:
        mem="15GB",
    threads: 4
    log:
        "logs/clustering/{subset}/plots/{assay}/{reduction_use}/correlation.log",


rule do_dot_plot:
    input:
        WorkflowResults(
            config["plotting"]["dot_plot"],
            "results/de/{subset}/plots/{assay}/dot_plot_by_{group}.png",
        ).create_path_list(),
    output:
        touch(temp("results/plotting/dot_plot_done")),
    localrule: True


use rule do_dot_plot as do_umap_plot with:
    input:
        WorkflowResults(
            config["plotting"]["umap_plot"],
            "results/clustering/{subset}/plots/{assay}/{reduction}/{group}_split_{split}.png",
        ).create_path_list(),
    output:
        touch(temp("results/plotting/umap_plot_done")),
    localrule: True


use rule do_dot_plot as do_correlation_plot with:
    input:
        WorkflowResults(
            config["plotting"]["umap_plot"],
            "results/clustering/{subset}/plots/{assay}/{reduction}/correlation.png",
        ).create_path_list(),
    output:
        touch(temp("results/plotting/cor_plot_done")),
    localrule: True


use rule do_dot_plot as do_de_plot with:
    input:
        expand(
            WorkflowResults(
                config["differential_expression"]["params"],
                "results/de/{subset}/plots/{assay}/{{plot_type}}_by_{group}.png",
            ).create_path_list(),
            plot_type=["heatmap", "de_plot"],
        ),
    output:
        touch(temp("results/plotting/de_plot_done")),
    localrule: True


use rule do_dot_plot as do_qc_plot with:
    input:
        expand("results/qc/violin/{type}_qc_plot.png", type=["pre", "post"]),
    output:
        touch(temp("results/plotting/qc_plot_done")),
    localrule: True


use rule do_dot_plot as do_bar_plot with:
    input:
        WorkflowResults(
            config["plotting"]["bar_plot"],
            "results/clustering/{subset}/plots/{assay}/{reduction}/barplot/{group}_by_{split}.png",
        ).create_path_list(),
    output:
        touch(temp("results/plotting/bar_plot_done")),
    localrule: True