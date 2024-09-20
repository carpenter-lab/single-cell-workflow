import os
pep = config["pep"]

rule gliph_prep:
    input:
        tcr=pep.sample_table.tcr_path,
    params:
        sample_name=pep.sample_table.sample_name,
        filter_chains=False,
        patient_id=pep.sample_table.patient_id,
        condition=pep.sample_table.condition,
        tool="gliph"
    output:
        tsv="results/gliph/input_tcr.tsv",
    log: "logs/gliph_prep.log"
    conda:
        "../envs/seurat.yml"
    script:
        "../scripts/tcr_metacluster_prep.R"


rule download_gliph:
    localrule: True
    output:
        "resources/gliph",
    shell:
        """
        DOWNLOAD_ROOT=http://50.255.35.37:8080/downloads/irtools

        if [[ $(uname) == "Linux" ]]; then
            URL="$DOWNLOAD_ROOT.centos"
        elif [[ $(uname) == "Darwin" ]]; then
            URL="$DOWNLOAD_ROOT.osx"
        else
            exit 1
        fi

        curl -so {output} "$URL"
        chmod +x {output}
        """


rule download_gliph_reference:
    localrule: True
    output:
        ref=directory("resources/gliph_ref"),
    shell:
        """
        temp=$(mktemp -t "gliph.XXXXXXXX.zip")
        temp_dir=$(mktemp -d)
        curl -so $temp http://50.255.35.37:8080/downloads/human_v2.0.zip
        unzip -q $temp -d $temp_dir
        mkdir -p {output}
        mv $temp_dir/*/* {output}
        """


rule render_gliph_template:
    input:
        "workflow/templates/gliph.txt",
    params:
        input_tcr=rules.gliph_prep.output.tsv,
        ref=rules.download_gliph_reference.output.ref,
        pwd=os.path.abspath(os.curdir),
        hla_file="data/HLA-ALL-Carpenter-forGLIPH2-2Aug2024.txt"
    output:
        temp("results/gliph/config.txt"),
    group:
        "configure_run_gliph"
    template_engine:
        "jinja2"


rule run_gliph:
    input:
        executable=rules.download_gliph.output,
        input_tcr=rules.gliph_prep.output.tsv,
        ref=rules.download_gliph_reference.output.ref,
        config_file=temp("results/gliph/config.txt"),
    output:
        prefix=directory("results/gliph/results"),
    group:
        "configure_run_gliph"
    script:
        "../scripts/run_gliph.sh"
