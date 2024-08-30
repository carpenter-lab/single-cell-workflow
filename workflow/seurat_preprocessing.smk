rule import_seurat:
    input:
        gex=pep.sample_table.gex_path,
        tcr=pep.sample_table.tcr_path,
    params:
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        sample_name=pep.sample_table.sample_name,
    output:
        seurat="results/import/object.rds",
    log:
        "logs/import.log",
    script:
        "scripts/import.R"


rule qc_setup:
    input:
        seurat=rules.import_seurat.output.seurat,
    output:
        seurat="results/qc/prep.rds",
    log:
        "logs/qc/setup.log",
    script:
        "scripts/quality_control_setup.R"


rule filter:
    input:
        seurat=rules.qc_setup.output.seurat,
    output:
        seurat="results/seurat_objects/all_data.rds",
    log:
        "logs/filter.log",
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
    script:
        "scripts/pca_reduction.R"
