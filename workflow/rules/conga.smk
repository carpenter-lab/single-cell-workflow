#!/usr/bin/env bash
import pandas as pd

pep = config["pep"]
rule install_conga:
    output:
        conga=directory("resources/conga"),
        setup="resources/conga/scripts/setup_10x_for_conga.py",
        merge="resources/conga/scripts/merge_samples.py",
        run="resources/conga/scripts/run_conga.py"
    conda: "../envs/conga.yml"
    shell:
        """
        git clone https://github.com/phbradley/conga.git {output.conga}
        
        cd {output.conga}/tcrdist_cpp || exit
        
        make

        cd .. || exit

        pip install -e .
        """


pep.sample_table["clones"] = [f"results/conga/clones/{sn}.tsv" for sn in pep.sample_table.sample_name]

rule setup_conga:
    input:
        tcr=pep.sample_table.tcr_path,
        setup=rules.install_conga.output.setup,
    output:
        clones=pep.sample_table.clones
    conda: "../envs/conga.yml"
    shell:
        """
        tcr=({input.tcr})
        clones=({output.clones})
        
        for i in "${{!tcr[@]}}"; do
            python {input.setup} \
                --filtered_contig_annotations_csvfile "${{tcr[i]}}"\
                --output_clones_file "${{clones[i]}}" \
                --organism human \
                --no_kpca 
        done
        """

rule merge_prep:
    input:
        clones=pep.sample_table.clones,
        gex=pep.sample_table.gex_path
    output:
        tsv="results/conga/merge.txt"
    script: "../scripts/conga_merge_file_prep.py"


rule merge:
    input:
        tsv=rules.merge_prep.output.tsv,
        merge=rules.install_conga.output.merge
    output:
        clones="results/conga/clones.tsv",
        gex="results/conga/gex.h5ad"
    conda: "../envs/conga.yml"
    shell:
        """
        python {input.merge} \
            --samples {input.tsv} \
            --output_clones_file {output.clones} \
            --output_gex_data {output.gex} \
            --organism human
        """

rule run_conga:
    input:
        gex=rules.merge.output.gex,
        clones=rules.merge.output.clones,
        exe=rules.install_conga.output.run
    output:
        x=directory("results/conga/results"),
        log="results/conga/results/x_log.txt"
    conda: "../envs/conga.yml"
    shell:
        """
        touch {output.log}
        python {input.exe} \
            --gex_data {input.gex} \
            --gex_data_type h5ad \
            --clones_file {input.clones}\
            --organism human \
            --all \
            --outfile_prefix {output.x}/x
        """