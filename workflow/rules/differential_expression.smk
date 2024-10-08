from common import (
    get_proper_clustering_output,
    validate_de_method,
    get_category_name,
    report_plot_labels,
    WorkflowResults,
)


rule differential_expression:
    input:
        seurat=get_proper_clustering_output(config, rules),
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        assay="{assay}",
        idents="{grouping}",
        test_use=validate_de_method(config["differential_expression"]["test"]),
    output:
        markers="results/de/{subset}/{assay}/{grouping}_markers.tsv",
    threads: lambda wildcards: 1 if wildcards.assay == "ADT" else 4
    resources:
        mem=lambda wildcards: "10G" if wildcards.assay == "ADT" else "40G",
        runtime="6h",
    log:
        "logs/de/{subset}/{assay}/{grouping}.log",
    conda:
        "../envs/differential_expression.yml"
    script:
        "../scripts/differential_expression.R"


rule table_to_html:
    input:
        table=rules.differential_expression.output.markers,
    output:
        html=report(
            "results/de/{subset}/{assay}/{grouping}_markers.html",
            category=get_category_name,
            subcategory="Differential Expression",
            labels=report_plot_labels,
        ),
        html_library=directory(
            "results/de/{subset}/{assay}/html_libs/{grouping}_markers"
        ),
    log:
        "logs/de/{subset}/{assay}/html_{grouping}.log",
    conda:
        "../envs/datatables_js.yml"
    script:
        "../scripts/dataframe_to_DataTables.R"


rule do_differential_expression:
    input:
        WorkflowResults(
            config["differential_expression"]["params"],
            "results/de/{subset}/{assay}/{group}_markers.{ext}",
            ext=["tsv", "html"],
        ).create_path_list(),
    output:
        touch(temp("results/de/done")),
    localrule: True
