from common import (
    get_proper_clustering_output,
    validate_de_method,
    get_category_name,
    report_plot_labels,
)


rule differential_expression:
    input:
        seurat=get_proper_clustering_output(config, rules),
    params:
        assay="{assay}",
        idents="{grouping}",
        test_use=validate_de_method(config["differential_expression"]["test"]),
    output:
        markers="results/de/{subset}/{assay}/{grouping}_markers.tsv",
    threads: lambda wildcards: 1 if wildcards.assay == "ADT" else min(workflow.cores / 2, 4)
    log:
        "logs/de/{subset}/{assay}/{grouping}.log",
    script:
        "scripts/differential_expression.R"


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
    script:
        "scripts/dataframe_to_DataTables.R"
