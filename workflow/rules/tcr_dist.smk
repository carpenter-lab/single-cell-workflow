pep = config["pep"]

rule tcrdist_prep:
    input:
        tcr=pep.sample_table.tcr_path,
    params:
        sample_name=pep.sample_table.sample_name,
        filter_chains=False,
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        tool="tcrdist"
    output:
        csv="results/tcr_dist/input_tcr.csv",
    log: "logs/tcrdist_prep.log"
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/tcr_metacluster_prep.R"

rule run_tcrdist:
    input:
        csv=rules.tcrdist_prep.output.csv,
    output:
        quasi_public="results/tcr_dist/quasi_public.csv",
        clone_df="results/tcr_dist/clone_df.csv",
        nn_summary="results/tcr_dist/nn_summary.csv",
    conda:
        "../envs/tcr_dist.yml"
    script:
        "../scripts/run_tcrdist.py"
