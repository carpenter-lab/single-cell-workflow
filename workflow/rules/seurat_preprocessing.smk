pep = config["pep"]


rule run_preprocessing:
    input:
        expand("results/qc/{type}_qc_plot.png", type=["pre", "post"]),
        expand(
            "results/pca/{subset}/object.rds",
            subset=config["subcluster"].get("all_data_key"),
        ),


rule import_rna:
    input:
        h5=pep.sample_table.gex_path,
    params:
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        sample_name=pep.sample_table.sample_name,
        h5_assay_name="Gene Expression",
    output:
        assay="results/import/rna.rds",
    log:
        "logs/import.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/import.R"


rule import_adt:
    input:
        h5=pep.sample_table.gex_path,
    params:
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        sample_name=pep.sample_table.sample_name,
        h5_assay_name="Antibody Capture",
    output:
        assay="results/import/adt.rds",
    log:
        "logs/import.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/import.R"


rule import_tcr:
    input:
        tcr=pep.sample_table.tcr_path,
    params:
        sample_name=pep.sample_table.sample_name,
        filter_chains=False,
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
    output:
        assay="results/import/tcr.rds",
    log:
        "logs/import_tcr.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/import_tcr.R"


def get_input_files(config):
    assays = config["preprocessing"]["assays"]
    assays = [assay.lower() for assay in assays]
    files = {f"{assay.upper()}": f"results/import/{assay}.rds" for assay in assays}
    return files


rule assemble_object:
    input:
        **get_input_files(config),
    output:
        seurat="results/import/object.rds",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/assemble_object.R"


rule qc_setup:
    input:
        seurat=rules.assemble_object.output,
    output:
        seurat="results/qc/prep.rds",
    log:
        "logs/qc/setup.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/quality_control_setup.R"


rule filter:
    input:
        seurat=rules.qc_setup.output.seurat,
    output:
        seurat="results/seurat_objects/all_data.rds",
    log:
        "logs/filter.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/filter.R"


rule normalise:
    input:
        seurat="results/seurat_objects/{subset}.rds",
    output:
        seurat="results/normalisation/{subset}.rds",
    params:
        method=config["preprocessing"]["normalisation"].get("method"),
        vars_to_regress=config["preprocessing"]["normalisation"]["vars_to_regress"],
    threads: 4
    resources:
        mem="18GB",
    log:
        "logs/normalise/{subset}.log",
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/normalise.R"


rule pca:
    input:
        seurat=rules.normalise.output.seurat,
    output:
        seurat="results/pca/{subset}/object.rds",
    threads: 4
    log:
        "logs/pca_{subset}.log",
    conda:
        (
            "../envs/trex.yml"
            if "TCR" in config["preprocessing"]["assays"]
            else "../envs/seurat.yml"
        )
    script:
        "../scripts/pca_reduction.R"
