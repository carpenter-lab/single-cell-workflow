#!/usr/bin/env bash

set -euo pipefail

# shellcheck disable=SC2154
EXEC="$PWD/${snakemake_input[executable]}"
CONFIG_FILE="$PWD/${snakemake_input[config_file]}"

# shellcheck disable=SC2154
mkdir -p "${snakemake_output[prefix]}"

cd "${snakemake_output[prefix]}" || exit

$EXEC -c "$CONFIG_FILE"

find x_* | sed "p;s/x_//g" | xargs -n2 mv
