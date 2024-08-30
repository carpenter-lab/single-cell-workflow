from common import get_subcluster_params


rule integrate:
    input:
        seurat=rules.pca.output.seurat,
    output:
        seurat="results/integration/{subset}.rds",
    threads: 4
    log:
        "logs/integration/{subset}.log",
    script:
        "scripts/integrate.R"


rule cluster:
    input:
        seurat=rules.integrate.output.seurat,
    params:
        assay="SCT",
        reductions=["pca", "harmony"],
    output:
        seurat="results/clustering/{subset}.rds",
    log:
        "logs/clustering/{subset}.log",
    script:
        "scripts/cluster.R"


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
    script:
        "scripts/azimuth.R"


rule subcluster:
    input:
        seurat=rules.run_azimuth.output.seurat,
    params:
        *get_subcluster_params(config),
    output:
        seurat="results/seurat_objects/{name}.rds",
    log:
        "logs/subclustering/{name}.log",
    script:
        "scripts/subcluster.R"
