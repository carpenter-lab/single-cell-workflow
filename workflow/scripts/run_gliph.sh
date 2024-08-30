#!/usr/bin/env bash

set -euo pipefail

# shellcheck disable=SC2154
touch "${snakemake_output[config_file]}"

cat <<- EOF > "${snakemake_output[config_file]}"
	out_prefix=x
	cdr3_file=$PWD/${snakemake_input[input_tcr]}
	refer_file=$PWD/${snakemake_input[ref]}/ref_CD4_v2.0.txt
	v_usage_freq_file=$PWD/${snakemake_input[ref]}/ref_V_CD4_v2.0.txt
	cdr3_length_freq_file=$PWD/${snakemake_input[ref]}/ref_L_CD4_v2.0.txt
	local_min_pvalue=0.001
	p_depth=1000
	global_convergence_cutoff=1
	simulation_depth=1000
	kmer_min_depth=3
	local_min_OVE=10
	algorithm=GLIPH2
	all_aa_interchangeable=1
EOF

EXEC="$PWD/${snakemake_input[executable]}"
CONFIG_FILE="$PWD/${snakemake_output[config_file]}"

mkdir -p "${snakemake_output[prefix]}"

cd "${snakemake_output[prefix]}" || exit

$EXEC -c "$CONFIG_FILE"

find x_* | sed "p;s/x_//g" | xargs -n2 mv
