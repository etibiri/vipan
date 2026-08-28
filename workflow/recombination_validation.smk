"""Snakemake rule for targeted RDP/3SEQ event validation.

Add this line to the project Snakefile:

    include: "workflow/recombination_validation.smk"
"""


VALIDATION_CONFIG = config.get("recombination_validation", {})


rule validate_recombination_events:
    input:
        events=VALIDATION_CONFIG.get(
            "events",
            "results/recombination/concordance_3seqraw_rdp_327.event_pairs.tsv",
        ),
        alignment=VALIDATION_CONFIG.get(
            "alignment",
            "results/alignment/splcv_aligned.fasta",
        ),
        script="scripts/validate_recombination_events.py",
    output:
        manifest="results/recombination/validation/validation_manifest.tsv",
        tree_summary="results/recombination/validation/segment_tree_summary.tsv",
        validation_calls="results/recombination/validation/event_validation_calls.tsv",
        done=touch("results/recombination/validation/validation.done"),
    log:
        "results/logs/validate_recombination_events.log",
    benchmark:
        "results/benchmarks/validate_recombination_events.tsv",
    conda:
        "../envs/recombination_validation.yaml"
    threads:
        VALIDATION_CONFIG.get("threads", 8)
    resources:
        mem_mb=VALIDATION_CONFIG.get("mem_mb", 16000),
        runtime=VALIDATION_CONFIG.get("runtime_minutes", 720),
    params:
        outdir="results/recombination/validation",
        tiers=VALIDATION_CONFIG.get("tiers", "Tier_B"),
        topology=VALIDATION_CONFIG.get("topology", "circular"),
        window=VALIDATION_CONFIG.get("window", 200),
        step=VALIDATION_CONFIG.get("step", 20),
        bootstraps=VALIDATION_CONFIG.get("bootstrap_replicates", 500),
        context_size=VALIDATION_CONFIG.get("context_size", 8),
        context_flag=(
            "--context-accessions " + VALIDATION_CONFIG["context_accessions"]
            if VALIDATION_CONFIG.get("context_accessions")
            else ""
        ),
        rdp_candidates=VALIDATION_CONFIG.get("rdp_candidates_per_role", 2),
        iqtree_flag=(
            "--run-iqtree"
            if VALIDATION_CONFIG.get("run_iqtree", True)
            else ""
        ),
        reuse_flag=(
            "--reuse-existing"
            if VALIDATION_CONFIG.get("reuse_existing", False)
            else ""
        ),
        formats=VALIDATION_CONFIG.get("formats", "png,pdf"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" results/logs results/benchmarks
        export MPLCONFIGDIR="${{TMPDIR:-/tmp}}/matplotlib-vipan"

        python "{input.script}" \
            --events "{input.events}" \
            --alignment "{input.alignment}" \
            --outdir "{params.outdir}" \
            --tiers "{params.tiers}" \
            --topology "{params.topology}" \
            --window {params.window} \
            --step {params.step} \
            --bootstrap-replicates {params.bootstraps} \
            --processes {threads} \
            --rdp-candidates-per-role {params.rdp_candidates} \
            --context-size {params.context_size} \
            {params.context_flag} \
            --formats "{params.formats}" \
            --iqtree-threads {threads} \
            {params.iqtree_flag} \
            {params.reuse_flag} \
            > "{log}" 2>&1
        """
