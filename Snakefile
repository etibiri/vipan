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
        os.path.join("data", "processed", "splcv_final.fasta"),
        os.path.join("data", "processed", "splcv_normalized.fasta"),
        os.path.join("results", "alignment", "splcv_aligned.fasta"),
        os.path.join("results", "analysis", "diversity_pi.tsv"),
        os.path.join("results", "plots", "pangenome_diversity.pdf")

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

# --- Rotate circular rule ---
rule normalize_genomes:
    """
    Step 4: Normalize circular genomes.
    1. Checks for conserved nonanucleotide (TAATATTAC).
    2. Flips to Reverse Complement if necessary.
    3. Rotates sequence to start at the motif.
    Ensures absolute co-linearity for Multiple Sequence Alignment.
    """
    input:
        fasta = "data/processed/splcv_final.fasta"
    output:
        fasta = "data/processed/splcv_normalized.fasta"
    params:
        motif = "TAATATTAC"
    conda:
        "envs/biopython.yaml"
    log:
        "results/logs/normalize_genomes.log"
    shell:
        """
        python scripts/rotate_sequences.py \
            --input {input.fasta} \
            --output {output.fasta} \
            --motif {params.motif} > {log} 2>&1
        """

################################################################################
# Multiple Sequence Alignment Step
################################################################################

rule align_genomes:
    """
    Step 5: Perform Multiple Sequence Alignment (MSA) using MAFFT.
    Running in a dedicated Conda environment.
    """
    input:
        fasta = "data/processed/splcv_normalized.fasta"
    output:
        alignment = "results/alignment/splcv_aligned.fasta"
    threads: 8 
    conda:
        "envs/mafft.yaml" 
    log:
        "results/logs/align_genomes.log"
    shell:
        """
        mkdir -p results/alignment
        mafft --auto --thread {threads} {input.fasta} > {output.alignment} 2> {log}
        """
rule analyze_diversity:
    """
    Step 6: Generate consensus sequence and calculate sliding-window diversity (Pi).
    Identifies variable regions and conserved motifs.
    """
    input:
        alignment = "results/alignment/splcv_aligned.fasta"
    output:
        consensus = "results/analysis/consensus.fasta",
        plot_data = "results/analysis/diversity_pi.tsv"
    conda:
        "envs/biopython.yaml"
    log:
        "results/logs/analyze_diversity.log"
    shell:
        """
        mkdir -p results/analysis
        python scripts/calculate_diversity.py \
            --input {input.alignment} \
            --consensus {output.consensus} \
            --plot_data {output.plot_data} > {log} 2>&1
        """

rule plot_diversity:
    """
    Step 7: Generate a PDF plot of Nucleotide Diversity (π).
    Visualizes hotspots of variation across the viral pangenome.
    """
    input:
        data = "results/analysis/diversity_pi.tsv"
    output:
        pdf = "results/plots/pangenome_diversity.pdf"
    conda:
        "envs/biopython.yaml"
    log:
        "results/logs/plot_diversity.log"
    shell:
        """
        mkdir -p results/plots
        python scripts/plot_diversity.py \
            --input {input.data} \
            --output {output.pdf} > {log} 2>&1
        """
