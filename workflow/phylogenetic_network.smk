"""Snakemake rules for phylogenetic-network construction (circular split
networks in the Neighbor-Net family; see scripts/build_phylogenetic_network.py
for the full method description and its explicit justification/limitations).

Builds one network per dataset already produced upstream in this workflow:
  - "329_dataset": results/<run_id>/alignment/initial_alignment.fasta
    (QC-passed genomes before the p-distance filter -- i.e. the sequence set
    referred to as "the 329-sequence dataset" in the manuscript)
  - "316_dataset": results/<run_id>/alignment/final_alignment.fasta
    (the p-distance-filtered set used for the primary RDP4/3SEQ analysis --
    "the 316-sequence dataset")

Pipeline (3 rules, in dependency order):
  1. distance_matrix       fasta -> pairwise p-distance TSV/PHYLIP only (fast;
                            exists purely to break the DAG cycle between step
                            2, which needs a distance matrix per dataset, and
                            step 3, which needs step 2's output table)
  2. classify_accessions    metadata + both distance matrices + existing
                            RDP4/3SEQ results -> accession_classification.tsv
  3. phylogenetic_network    fasta + accession_classification.tsv -> full
                            circular split network (recomputes its own
                            distance matrix internally; step 1's output is
                            not reused here, it is only an input to step 2)

Add this line to the project Snakefile to wire these targets into the
default DAG (matching the existing opt-in pattern used by
workflow/recombination_validation.smk and workflow/filtering_sensitivity.smk):

    include: "workflow/phylogenetic_network.smk"

Runnable standalone, e.g. on the existing SLURM profile:

    snakemake --profile profiles/slurm phylogenetic_networks_all

or locally:

    snakemake -c4 --use-conda phylogenetic_networks_all
"""

NETWORK_CONFIG = config.get("phylogenetic_network", {})

NETWORK_DATASETS = {
    "329_dataset": {
        "fasta": f"{RUN_DIR}/alignment/initial_alignment.fasta",
        "prefix": "splcv_329",
        "title": "Sweepovirus circular split network -- full QC-passed dataset (n=329)",
    },
    "316_dataset": {
        "fasta": f"{RUN_DIR}/alignment/final_alignment.fasta",
        "prefix": "splcv_316",
        "title": "Sweepovirus circular split network -- distance-filtered dataset (n=316)",
    },
}


rule distance_matrix:
    """Fast, network-plotting-independent pairwise p-distance matrix for one
    dataset (used only as an input to classify_accessions, see module
    docstring)."""
    input:
        fasta=lambda wc: NETWORK_DATASETS[wc.dataset]["fasta"],
        script="scripts/build_phylogenetic_network.py",
    output:
        pairwise=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.pairwise_distances.tsv",
        phylip=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.distance_matrix.phylip",
    log:
        f"{RUN_DIR}/logs/distance_matrix_{{dataset}}_{{prefix}}.log",
    conda:
        "../envs/recombination_core.yaml"
    params:
        outdir=lambda wc: f"{RUN_DIR}/phylogenetic_network/{wc.dataset}",
        prefix=lambda wc: NETWORK_DATASETS[wc.dataset]["prefix"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        python "{input.script}" \
            --fasta "{input.fasta}" \
            --outdir "{params.outdir}" \
            --prefix "{params.prefix}" \
            --distances-only \
            > "{log}" 2>&1
        """


rule classify_accessions:
    """Build the accession -> species / phylogenetic-cluster / genetic-group /
    recombinant-status classification table used to annotate both networks.
    Terminology is restricted to species / genetic group / lineage / cluster;
    'subspecies' is never used (not a recognized rank for this genus)."""
    input:
        metadata=config["input"].get("metadata", "data/raw/splcv_all-0726.tsv"),
        dist_329="{run}/phylogenetic_network/329_dataset/{prefix}.pairwise_distances.tsv".format(
            run=RUN_DIR, prefix=NETWORK_DATASETS["329_dataset"]["prefix"]
        ),
        dist_316="{run}/phylogenetic_network/316_dataset/{prefix}.pairwise_distances.tsv".format(
            run=RUN_DIR, prefix=NETWORK_DATASETS["316_dataset"]["prefix"]
        ),
        rdp_robust=f"{RUN_DIR}/recombination/rdp_events_robust.tsv",
        strict_calls=f"{RUN_DIR}/validation/event_validation_calls_strict.tsv",
        script="scripts/classify_accessions.py",
    output:
        classification=f"{RUN_DIR}/phylogenetic_network/accession_classification.tsv",
    log:
        f"{RUN_DIR}/logs/classify_accessions.log",
    conda:
        "../envs/recombination_core.yaml"
    params:
        coarse_threshold=NETWORK_CONFIG.get("coarse_cluster_p_distance", 0.09),
        split_size=NETWORK_CONFIG.get("lineage_split_size", 25),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.classification}")" "$(dirname "{log}")"
        python "{input.script}" \
            --metadata "{input.metadata}" \
            --distances-329 "{input.dist_329}" \
            --distances-316 "{input.dist_316}" \
            --rdp-robust "{input.rdp_robust}" \
            --strict-calls "{input.strict_calls}" \
            --coarse-threshold {params.coarse_threshold} \
            --split-size {params.split_size} \
            --output "{output.classification}" \
            > "{log}" 2>&1
        """


rule phylogenetic_network:
    """Full circular split network for one dataset (329_dataset or
    316_dataset), coloured by the classification table's lineage column."""
    input:
        fasta=lambda wc: NETWORK_DATASETS[wc.dataset]["fasta"],
        groups=f"{RUN_DIR}/phylogenetic_network/accession_classification.tsv",
        script="scripts/build_phylogenetic_network.py",
    output:
        summary=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.network_summary.json",
        splits=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.splits.tsv",
        nexus=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.splits.nexus",
        order=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.circular_order.tsv",
        pdf=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.network.pdf",
        svg=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.network.svg",
        png=f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.network.png",
    log:
        f"{RUN_DIR}/logs/phylogenetic_network_{{dataset}}_{{prefix}}.log",
    benchmark:
        f"{RUN_DIR}/benchmarks/phylogenetic_network_{{dataset}}_{{prefix}}.tsv",
    conda:
        "../envs/recombination_core.yaml"
    threads: 1
    resources:
        mem_mb=NETWORK_CONFIG.get("mem_mb", 8000),
        runtime=NETWORK_CONFIG.get("runtime_minutes", 60),
    params:
        outdir=lambda wc: f"{RUN_DIR}/phylogenetic_network/{wc.dataset}",
        prefix=lambda wc: NETWORK_DATASETS[wc.dataset]["prefix"],
        title=lambda wc: NETWORK_DATASETS[wc.dataset]["title"],
        max_split_side=NETWORK_CONFIG.get("max_split_side", 5),
        min_weight_fraction=NETWORK_CONFIG.get("min_split_weight_fraction", 0.002),
        group_column=lambda wc: f"lineage_{wc.dataset.split('_')[0]}",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        python "{input.script}" \
            --fasta "{input.fasta}" \
            --outdir "{params.outdir}" \
            --prefix "{params.prefix}" \
            --title "{params.title}" \
            --max-split-side {params.max_split_side} \
            --min-split-weight-fraction {params.min_weight_fraction} \
            --groups "{input.groups}" \
            --group-column "{params.group_column}" \
            > "{log}" 2>&1
        """


rule phylogenetic_networks_all:
    input:
        expand(
            f"{RUN_DIR}/phylogenetic_network/{{dataset}}/{{prefix}}.network_summary.json",
            zip,
            dataset=list(NETWORK_DATASETS.keys()),
            prefix=[v["prefix"] for v in NETWORK_DATASETS.values()],
        ),
