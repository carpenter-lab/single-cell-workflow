
rule prep:
    input:
        seurat="results/seurat_objects/all_data.rds",
    output:
        csv="results/tcr_dist/input_tcr.csv",
    conda: "../envs/seurat.yml"
    script:
        "../scripts/tcrdist_prep.R"

rule run_tcrdist:
    input:
        csv=rules.prep.output.csv
    output:
        quasi_public="results/tcr_dist/quasi_public.csv",
        clone_df="results/tcr_dist/clone_df.csv",
        nn_summary="results/tcr_dist/nn_summary.csv"
    conda: "../envs/tcr_dist.yml"
    script:
        "../scripts/run_tcrdist.py"