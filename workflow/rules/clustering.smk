from typing import Callable

from common import get_subcluster_params, valid_dict_key


rule integrate:
    input:
        seurat="results/pca/{subset}/object.rds",
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        method="harmony",
    output:
        seurat="results/integration/{subset}.rds",
    threads: 4
    log:
        "logs/integration/{subset}.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/integrate.R"


rule cluster:
    input:
        seurat="results/integration/{subset}.rds",
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        assay="SCT",
        reductions=["pca", "harmony"],
    output:
        seurat="results/clustering/{subset}.rds",
    log:
        "logs/clustering/{subset}.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/cluster.R"


rule label_clusters:
    input:
        seurat=expand(
            rules.cluster.output.seurat,
            subset=config["cluster"].get("all_data_key"),
        ),
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        **config["cluster"].get("labels"),
    output:
        seurat=f"results/clustering/all_data/{config['cluster'].get('labels').get('new_group_by')}.rds",
    conda:
        "../envs/seurat.yml"
    log:
        "logs/clustering/label_clusters.log",
    script:
        "../scripts/label_clusters.R"


def decide_label_input(config: dict) -> Callable:
    def _decide_label_input(wildcards: dict) -> str:
        if valid_dict_key(config["cluster"], "labels"):
            if not config["cluster"]["labels"]["skip"]:
                return f"results/clustering/all_data/{config['cluster'].get('labels').get('new_group_by')}.rds"
        return "results/clustering/{subset}.rds".format(subset=config["cluster"].get("all_data_key"))

    return _decide_label_input


rule annotate_tcr_seurat:
    input:
        seurat=decide_label_input(config),
        matrix_dir=(
            "results/import/bpcells_backing"
            if config["preprocessing"]["use_bpcells"]
            else None
        ),
    params:
        tcr_patterns={
            "TB": [
                "SLG%WET",
                "SSPGQQGG%NYG",
                "S%GTESNQP",
                "SPGR%E",
                "SYSGR%TE",
                "SLG%LE",
                "SQG%AYNE",
                "SRG%QP",
                "G%GEGQP",
                "SQER%YG",
                "SQEGR%NQP",
                "SLGTESN%P",
                "SPGTSG%DT",
                "S%VTSGTYE",
                "RTG%YE",
            ],
            "Other": "SRDR%SYG",
        },
    output:
        seurat="results/gliph/labelled_tcr.rds",
    log:
        "logs/gliph/label_tcr.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/label_tcr.R"


if valid_dict_key(config["cluster"], "subclusters"):

    rule subcluster:
        input:
            seurat=rules.annotate_tcr_seurat.output.seurat,
            matrix_dir=(
                "results/import/bpcells_backing"
                if config["preprocessing"]["use_bpcells"]
                else None
            ),
        params:
            *get_subcluster_params(config),
        output:
            seurat="results/seurat_objects/{name}.rds",
        log:
            "logs/subclustering/{name}.log",
        conda:
            "../envs/seurat.yml"
        script:
            "../scripts/subcluster.R"
