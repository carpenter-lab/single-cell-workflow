rule import_seurat:
    input:
        gex=pep.sample_table.gex_path,
    params:
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        sample_name=pep.sample_table.sample_name,
    output:
        seurat="results/import/object.rds",
    log:
        "logs/import.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/import.R"

rule import_tcr:
    input:
        seurat=rules.import_seurat.output.seurat,
        tcr=pep.sample_table.tcr_path,
    params:
        sample_name=pep.sample_table.sample_name,
        filter_chains=False
    output:
        seurat="results/import/object_with_tcr.rds",
    log:
        "logs/import_tcr.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/import_tcr.R"

rule qc_setup:
    input:
        seurat=rules.import_tcr.output.seurat,
    output:
        seurat="results/qc/prep.rds",
    log:
        "logs/qc/setup.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/quality_control_setup.R"


rule filter:
    input:
        seurat=rules.qc_setup.output.seurat,
    output:
        seurat="results/seurat_objects/all_data.rds",
    log:
        "logs/filter.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/filter.R"


rule normalise:
    input:
        seurat="results/seurat_objects/{subset}.rds",
    output:
        seurat="results/normalisation/{subset}.rds",
    params:
        method=config["preprocessing"]["normalisation"].get("method"),
        vars_to_regress=config["preprocessing"]["normalisation"]["vars_to_regress"],
    threads: 4
    log:
        "logs/normalise/{subset}.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/normalise.R"


rule pca:
    input:
        seurat=rules.normalise.output.seurat,
    output:
        seurat="results/pca/{subset}/object.rds",
    threads: 4
    log:
        "logs/pca_{subset}.log",
    conda: "envs/seurat.yml"
    script:
        "scripts/pca_reduction.R"
