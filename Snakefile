################################################################################
# Recombination Analysis — SWEEPOVIRUSES
# Author: Ezechiel B. TIBIRI
# Reproducible, offline-ready Snakemake workflow.
# Uses existing Conda envs (envs/) and provided scripts (scripts/).
################################################################################

import os

# Load configuration parameters
configfile: "config.yaml"

# --- Target Rule ---
# This rule tells Snakemake what the final desired output is
rule all:
    input:
        os.path.join("data", "raw", "splcv_all.fasta"),
        os.path.join("data", "raw", "splcv_metadata.tsv")

# --- Workflow Rules ---

rule download_genomes:
    """
    Step 1: Fetch sequences and metadata from NCBI Nucleotide.
    Uses the custom python script in scripts/ and the parameters from config.yaml
    """
    output:
        fasta = "data/raw/splcv_all.fasta",
        metadata = "data/raw/splcv_metadata.tsv"
    params:
        email = config["ncbi"]["email"],
        min_len = config["filters"]["min_len"],
        max_len = config["filters"]["max_len"]
    conda:
        "envs/biopython.yaml"
    log:
        "results/logs/download_genomes.log"
    shell:
        """
        python scripts/fetch_splcv.py \
            --email {params.email} \
            --min_len {params.min_len} \
            --max_len {params.max_len} \
            --fasta {output.fasta} \
            --tsv {output.metadata} > {log} 2>&1
        """
