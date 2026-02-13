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
        os.path.join("data", "raw", "splcv_metadata.tsv"),
        os.path.join("data", "processed", "splcv_clean.fasta"),
        os.path.join("data", "processed", "splcv_final.fasta")

# --- Workflow Rules ---

# --- Download sequences ---
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
        python scripts/fetch_sweepo.py \
            --email {params.email} \
            --min_len {params.min_len} \
            --max_len {params.max_len} \
            --fasta {output.fasta} \
            --tsv {output.metadata} > {log} 2>&1
        """

# --- N quality filter ---
rule filter_quality:
    """
    Step 2: Remove sequences with more than 1% of Ns to ensure pangenome quality.
    """
    input:
        fasta = "data/raw/splcv_all.fasta"
    output:
        fasta = "data/processed/splcv_clean.fasta"
    params:
        max_n = 1.0  # Threshold: 1% max N
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/filter_n.py \
            --input {input.fasta} \
            --output {output.fasta} \
            --max_n_pct {params.max_n}
        """
# --- Fix duplicates rule ---
rule remove_duplicates:
    """
    Step 3: Remove 100% identical sequences (dereplication).
    Reduces redundancy and prevents frequency bias in the pangenome.
    """
    input:
        fasta = "data/processed/splcv_clean.fasta"
    output:
        fasta = "data/processed/splcv_final.fasta"
    conda:
        "envs/biopython.yaml"
    log:
        "results/logs/remove_duplicates.log"
    shell:
        """
        python scripts/remove_duplicates.py \
            --input {input.fasta} \
            --output {output.fasta} > {log} 2>&1
        """
