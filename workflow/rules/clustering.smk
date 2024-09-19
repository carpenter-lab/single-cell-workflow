from common import get_subcluster_params


rule integrate:
    input:
        seurat="results/pca/{subset}/object.rds",
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


rule run_azimuth:
    input:
        seurat=expand(
            rules.cluster.output.seurat,
            subset=config["subcluster"].get("all_data_key"),
        ),
    params:
        reference="lungref",
    output:
        seurat="results/clustering/azimuth_annotation.rds",
    log:
        "logs/azimuth.log",
    threads: 4
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/azimuth.R"


rule annotate_tcr_seurat:
    input:
        seurat=expand(
            rules.cluster.output.seurat,
            subset=config["subcluster"].get("all_data_key"),
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


rule subcluster:
    input:
        seurat=rules.annotate_tcr_seurat.output.seurat,
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
