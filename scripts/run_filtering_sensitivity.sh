#!/usr/bin/env bash
set -euo pipefail

cores="${1:-8}"
script_directory=$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)
project_directory=$(
    cd "${script_directory}/.."
    pwd
)

cd "${project_directory}"

snakemake \
    --use-conda \
    --cores "${cores}" \
    --configfile configs/medoid_sensitivity.yaml \
    prepare_filtered_alignment

snakemake \
    --use-conda \
    --cores "${cores}" \
    --configfile configs/reference_NC_004650_sensitivity.yaml \
    prepare_filtered_alignment

snakemake \
    --snakefile workflow/filtering_sensitivity.smk \
    --use-conda \
    --cores 1
