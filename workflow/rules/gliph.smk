
rule gliph_prep1:
    input:
        seurat=rules.filter.output.seurat,
    output:
        tsv="results/gliph/input_tcr.tsv",
    conda: "envs/seurat.yml"
    script:
        "scripts/gliph_prep.R"


rule download_gliph:
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
    output:
        directory("resources/gliph_ref"),
    shell:
        """
        temp=$(mktemp -t "gliph.XXXXXXXX.zip")
        temp_dir=$(mktemp -d)
        curl -so $temp http://50.255.35.37:8080/downloads/human_v2.0.zip
        unzip -q $temp -d $temp_dir
        mkdir -p {output}
        mv $temp_dir/*/* {output}
        """


rule run_gliph:
    input:
        executable=rules.download_gliph.output,
        input_tcr=rules.gliph_prep1.output.tsv,
        ref=rules.download_gliph_reference.output,
    output:
        prefix=directory("results/gliph/results"),
        config_file=temp("results/gliph/config.txt"),
    script:
        "scripts/run_gliph.sh"
