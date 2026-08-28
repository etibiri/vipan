"""End-to-end, audited SPLCV recombination workflow.

The workflow has one deliberate external boundary: RDP4 GUI analysis.  The
target ``prepare_recombination_inputs`` creates the final alignment and the
standalone 3SEQ result.  Export the RDP CSV to the path configured under
``rdp.report_csv`` and then run the default ``all`` target.
"""

configfile: "config.yaml"

from pathlib import Path


RUN_ID = config["run_id"]
RUN_DIR = f"results/{RUN_ID}"
PRE = config["preprocessing"]
ALN = config["alignment"]
DIST = config["distance_filter"]
THREESEQ = config["threeseq"]
RDP = config["rdp"]
COMP = config["comparison"]
VALID = config["validation"]
SIGNAL = config["phylogenetic_signal"]


def bool_flag(condition, flag):
    return flag if condition else ""


rule all:
    input:
        f"{RUN_DIR}/recombination/concordance.summary.tsv",
        f"{RUN_DIR}/validation/event_validation_summary_strict.tsv",
        f"{RUN_DIR}/validation/phylogenetic_signal/event_phylogenetic_signal_summary.tsv",
        f"{RUN_DIR}/provenance/run_manifest.json",


rule prepare_recombination_inputs:
    input:
        f"{RUN_DIR}/alignment/final_alignment.fasta",
        f"{RUN_DIR}/qc/final_alignment_summary.tsv",
        f"{RUN_DIR}/recombination/3seq_output.rec",


rule prepare_filtered_alignment:
    input:
        f"{RUN_DIR}/alignment/final_alignment.fasta",
        f"{RUN_DIR}/qc/final_alignment_summary.tsv",
        f"{RUN_DIR}/qc/distance_filter_decisions.tsv",
        f"{RUN_DIR}/qc/distance_filter.pdf",


rule prepare_sequences:
    input:
        fasta=config["input"]["fasta"],
        script="scripts/prepare_sequences.py",
    output:
        fasta=f"{RUN_DIR}/sequences/prepared.fasta",
        excluded=f"{RUN_DIR}/sequences/pre_alignment_excluded.fasta",
        manifest=f"{RUN_DIR}/qc/sequence_preparation_manifest.tsv",
    log:
        f"{RUN_DIR}/logs/prepare_sequences.log",
    conda:
        "envs/recombination_core.yaml"
    params:
        min_length=PRE.get("min_length", 2400),
        max_ambiguous=PRE.get("max_ambiguous_fraction", 0.05),
        motif=PRE.get("circular_motif", "TAATATTAC"),
        motif_mismatches=PRE.get("max_motif_mismatches", 1),
        orient_flag=bool_flag(PRE.get("orient_by_motif", True), "--orient-by-motif"),
        require_motif_flag=bool_flag(PRE.get("require_motif", False), "--require-motif"),
        deduplicate_flag=bool_flag(PRE.get("deduplicate_exact", True), "--deduplicate-exact"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.fasta}")" "$(dirname "{output.manifest}")" "$(dirname "{log}")"
        python "{input.script}" \
            --input "{input.fasta}" \
            --output-fasta "{output.fasta}" \
            --excluded-fasta "{output.excluded}" \
            --manifest "{output.manifest}" \
            --min-length {params.min_length} \
            --max-ambiguous-fraction {params.max_ambiguous} \
            --circular-motif "{params.motif}" \
            --max-motif-mismatches {params.motif_mismatches} \
            {params.orient_flag} \
            {params.require_motif_flag} \
            {params.deduplicate_flag} \
            > "{log}" 2>&1
        """

rule initial_alignment:
    input:
        fasta=f"{RUN_DIR}/sequences/prepared.fasta",
    output:
        alignment=f"{RUN_DIR}/alignment/initial_alignment.fasta",
    log:
        f"{RUN_DIR}/logs/initial_alignment.log",
    benchmark:
        f"{RUN_DIR}/benchmarks/initial_alignment.tsv",
    threads:
        ALN.get("threads", 8)
    conda:
        "envs/recombination_core.yaml"
    params:
        options=ALN.get("mafft_options", "--auto"),
        benchmark_dir=f"{RUN_DIR}/benchmarks",
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname "{output.alignment}")" \
            "$(dirname "{log}")" \
            "{params.benchmark_dir}"

        test -s "{input.fasta}"

        mafft \
            {params.options} \
            --thread {threads} \
            "{input.fasta}" \
            > "{output.alignment}" \
            2> "{log}"

        test -s "{output.alignment}"
        """

rule filter_alignment_outliers:
    input:
        alignment=f"{RUN_DIR}/alignment/initial_alignment.fasta",
        script="scripts/filter_alignment_outliers.py",
    output:
        unaligned=f"{RUN_DIR}/sequences/distance_filtered.fasta",
        aligned=f"{RUN_DIR}/alignment/initial_alignment_kept.fasta",
        excluded=f"{RUN_DIR}/alignment/initial_alignment_excluded.fasta",
        metrics=f"{RUN_DIR}/qc/distance_filter_decisions.tsv",
        distances=f"{RUN_DIR}/qc/pairwise_distances_initial.tsv",
        summary=f"{RUN_DIR}/qc/distance_filter_summary.json",
        plot_png=f"{RUN_DIR}/qc/distance_filter.png",
        plot_pdf=f"{RUN_DIR}/qc/distance_filter.pdf",
    log:
        f"{RUN_DIR}/logs/filter_alignment_outliers.log",
    conda:
        "envs/recombination_core.yaml"
    params:
        plot_prefix=f"{RUN_DIR}/qc/distance_filter",
        mode=DIST.get("mode", "medoid"),
        reference_flag=(
            f"--reference-accession '{DIST['reference_accession']}'"
            if DIST.get("reference_accession")
            else ""
        ),
        max_distance=DIST.get("max_p_distance", 0.25),
        min_comparable=DIST.get("min_comparable_sites", 1000),
        min_overlap=DIST.get("min_overlap_fraction", 0.60),
        min_valid=DIST.get("min_valid_bases", 2400),
        max_gap=DIST.get("max_gap_fraction", 0.20),
        max_ambiguous=DIST.get("max_ambiguity_fraction", 0.05),
        keep_flag=(f"--manual-keep '{DIST['manual_keep']}'" if DIST.get("manual_keep") else ""),
        exclude_flag=(f"--manual-exclude '{DIST['manual_exclude']}'" if DIST.get("manual_exclude") else ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.unaligned}")" "$(dirname "{output.metrics}")" "$(dirname "{log}")"
        export MPLCONFIGDIR="${{TMPDIR:-/tmp}}/matplotlib-vipan-filter"
        python "{input.script}" \
            --alignment "{input.alignment}" \
            --kept-fasta "{output.unaligned}" \
            --kept-aligned-fasta "{output.aligned}" \
            --excluded-fasta "{output.excluded}" \
            --metrics "{output.metrics}" \
            --pairwise-distances "{output.distances}" \
            --summary-json "{output.summary}" \
            --plot-prefix "{params.plot_prefix}" \
            --mode "{params.mode}" \
            --max-p-distance {params.max_distance} \
            --min-comparable-sites {params.min_comparable} \
            --min-overlap-fraction {params.min_overlap} \
            --min-valid-bases {params.min_valid} \
            --max-gap-fraction {params.max_gap} \
            --max-ambiguity-fraction {params.max_ambiguous} \
            {params.reference_flag} \
            {params.keep_flag} \
            {params.exclude_flag} \
            > "{log}" 2>&1
        """

rule final_alignment:
    input:
        fasta=f"{RUN_DIR}/sequences/distance_filtered.fasta",
    output:
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
    log:
        f"{RUN_DIR}/logs/final_alignment.log",
    benchmark:
        f"{RUN_DIR}/benchmarks/final_alignment.tsv",
    threads:
        ALN.get("threads", 8)
    conda:
        "envs/recombination_core.yaml"
    params:
        options=ALN.get("mafft_options", "--auto"),
        benchmark_dir=f"{RUN_DIR}/benchmarks",
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname "{output.alignment}")" \
            "$(dirname "{log}")" \
            "{params.benchmark_dir}"

        test -s "{input.fasta}"

        mafft \
            {params.options} \
            --thread {threads} \
            "{input.fasta}" \
            > "{output.alignment}" \
            2> "{log}"

        test -s "{output.alignment}"
        """

rule summarize_final_alignment:
    input:
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        script="scripts/summarize_alignment.py",
    output:
        summary=f"{RUN_DIR}/qc/final_alignment_summary.tsv",
        per_sequence=f"{RUN_DIR}/qc/final_alignment_per_sequence.tsv",
        json=f"{RUN_DIR}/qc/final_alignment_summary.json",
    log:
        f"{RUN_DIR}/logs/summarize_final_alignment.log",
    conda:
        "envs/recombination_core.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.summary}")" "$(dirname "{log}")"
        python "{input.script}" \
            --alignment "{input.alignment}" \
            --summary-tsv "{output.summary}" \
            --per-sequence-tsv "{output.per_sequence}" \
            --summary-json "{output.json}" \
            > "{log}" 2>&1
        """


rule run_3seq:
    input:
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        executable=THREESEQ["executable"],
        ptable=THREESEQ["ptable"],
    output:
        rec=f"{RUN_DIR}/recombination/3seq_output.rec",
        native_log=f"{RUN_DIR}/recombination/3seq_native.3s.log",
    log:
        f"{RUN_DIR}/logs/run_3seq.log",
    benchmark:
        f"{RUN_DIR}/benchmarks/run_3seq.tsv",
    threads:
        THREESEQ.get("threads", 1)
    resources:
        mem_mb=THREESEQ.get("mem_mb", 2000),
        runtime=THREESEQ.get("runtime_minutes", 720),
    params:
        run_id=THREESEQ.get("run_id", "splcv_3seq"),
        benchmark_dir=f"{RUN_DIR}/benchmarks",
    shell:
        r"""
        set -euo pipefail

        project_dir=$(pwd)

        alignment=$(realpath "{input.alignment}")
        executable=$(realpath "{input.executable}")
        ptable=$(realpath "{input.ptable}")

        output_rec="$project_dir/{output.rec}"
        output_native_log="$project_dir/{output.native_log}"
        run_log="$project_dir/{log}"

        mkdir -p \
            "$(dirname "$output_rec")" \
            "$(dirname "$output_native_log")" \
            "$(dirname "$run_log")" \
            "{params.benchmark_dir}"

        if [ ! -s "$alignment" ]; then
            echo "ERROR: alignment absent ou vide : $alignment" >&2
            exit 1
        fi

        if [ ! -x "$executable" ]; then
            echo "ERROR: exécutable 3SEQ absent ou non exécutable : $executable" >&2
            exit 1
        fi

        if [ ! -s "$ptable" ]; then
            echo "ERROR: table de probabilités absente ou vide : $ptable" >&2
            exit 1
        fi

        workdir=$(mktemp -d "${{TMPDIR:-/tmp}}/3seq_run.XXXXXX")

        cleanup() {{
            rm -rf "$workdir"
        }}

        trap cleanup EXIT

        cd "$workdir"

        printf 'Y\nY\nY\n' |
            "$executable" \
                -full "$alignment" \
                -ptable "$ptable" \
                -id "{params.run_id}" \
                > "$run_log" 2>&1

        if grep -q "ERROR:" "$run_log"; then
            cat "$run_log" >&2
            exit 1
        fi

        if [ ! -s "{params.run_id}.3s.rec" ]; then
            echo \
                "ERROR: 3SEQ n'a pas produit le fichier {params.run_id}.3s.rec" \
                >&2
            cat "$run_log" >&2
            exit 1
        fi

        if [ ! -s "{params.run_id}.3s.log" ]; then
            echo \
                "ERROR: 3SEQ n'a pas produit le fichier {params.run_id}.3s.log" \
                >&2
            cat "$run_log" >&2
            exit 1
        fi

        mv \
            "{params.run_id}.3s.rec" \
            "$output_rec"

        mv \
            "{params.run_id}.3s.log" \
            "$output_native_log"

        test -s "$output_rec"
        test -s "$output_native_log"
        """

rule parse_rdp_report:
    input:
        report=RDP["report_csv"],
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        script="scripts/parse_rdp_csv.py",
    output:
        all_events=f"{RUN_DIR}/recombination/rdp_events_all.tsv",
        robust=f"{RUN_DIR}/recombination/rdp_events_robust.tsv",
    log:
        f"{RUN_DIR}/logs/parse_rdp_report.log",
    conda:
        "envs/recombination_core.yaml"
    params:
        min_methods=RDP.get("min_methods", 3),
        min_tested=RDP.get("min_tested_methods", 7),
        min_support=RDP.get("min_support_fraction", 0.50),
        max_missing=RDP.get("max_missing_methods", 2),
        p_cutoff=RDP.get("method_pvalue_cutoff", 0.05),
        min_fragment=RDP.get("min_fragment_length", 100),
        internal_flag=bool_flag(RDP.get("include_internal_3seq", False), "--include-internal-3seq"),
        known_flag=bool_flag(RDP.get("require_known_parents", True), "--require-known-parents"),
        caret_flag=bool_flag(RDP.get("allow_caret", False), "--allow-caret"),
        tilde_flag=bool_flag(RDP.get("allow_tilde", False), "--allow-tilde"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.all_events}")" "$(dirname "{log}")"
        python "{input.script}" \
            --input "{input.report}" \
            --alignment "{input.alignment}" \
            --output-all "{output.all_events}" \
            --output-robust "{output.robust}" \
            --min-methods {params.min_methods} \
            --min-tested-methods {params.min_tested} \
            --min-support-fraction {params.min_support} \
            --max-missing-methods {params.max_missing} \
            --method-pvalue-cutoff {params.p_cutoff} \
            --min-fragment-length {params.min_fragment} \
            {params.internal_flag} \
            {params.known_flag} \
            {params.caret_flag} \
            {params.tilde_flag} \
            > "{log}" 2>&1
        """


rule compare_3seq_rdp_events:
    input:
        threeseq=f"{RUN_DIR}/recombination/3seq_output.rec",
        rdp=f"{RUN_DIR}/recombination/rdp_events_robust.tsv",
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        script="scripts/compare_3seq_rdp_events.py",
    output:
        pairs=f"{RUN_DIR}/recombination/concordance.event_pairs.tsv",
        recombinants=f"{RUN_DIR}/recombination/concordance.recombinant_summary.tsv",
        summary=f"{RUN_DIR}/recombination/concordance.summary.tsv",
    log:
        f"{RUN_DIR}/logs/compare_3seq_rdp_events.log",
    conda:
        "envs/recombination_core.yaml"
    params:
        prefix=f"{RUN_DIR}/recombination/concordance",
        topology=COMP.get("topology", "circular"),
        coordinate_system=COMP.get("rdp_coordinate_system", "alignment"),
        tolerance=COMP.get("tolerance", 50),
        multiplier=COMP.get("moderate_multiplier", 2),
        parent_penalty=COMP.get("parent_penalty", 0),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.pairs}")" "$(dirname "{log}")"
        python "{input.script}" \
            --threeseq "{input.threeseq}" \
            --rdp "{input.rdp}" \
            --alignment "{input.alignment}" \
            --topology "{params.topology}" \
            --rdp-coordinate-system "{params.coordinate_system}" \
            --tolerance {params.tolerance} \
            --moderate-multiplier {params.multiplier} \
            --parent-penalty {params.parent_penalty} \
            --out-prefix "{params.prefix}" \
            > "{log}" 2>&1
        """


rule validate_recombination_events:
    input:
        events=f"{RUN_DIR}/recombination/concordance.event_pairs.tsv",
        alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        script="scripts/validate_recombination_events.py",
    output:
        manifest=f"{RUN_DIR}/validation/validation_manifest.tsv",
        tree_summary=f"{RUN_DIR}/validation/segment_tree_summary.tsv",
        calls=f"{RUN_DIR}/validation/event_validation_calls.tsv",
    log:
        f"{RUN_DIR}/logs/validate_recombination_events.log",
    benchmark:
        f"{RUN_DIR}/benchmarks/validate_recombination_events.tsv",
    threads:
        VALID.get("threads", 8)
    resources:
        mem_mb=VALID.get("mem_mb", 16000),
        runtime=VALID.get("runtime_minutes", 720),
    conda:
        "envs/recombination_validation.yaml"
    params:
        outdir=f"{RUN_DIR}/validation",
        benchmark_dir=f"{RUN_DIR}/benchmarks",
        tiers=VALID.get("tiers", "Tier_B"),
        topology=COMP.get("topology", "circular"),
        window=VALID.get("window", 200),
        step=VALID.get("step", 20),
        bootstrap=VALID.get("bootstrap_replicates", 500),
        seed=VALID.get("seed", 20260820),
        rdp_candidates=VALID.get(
            "rdp_candidates_per_role",
            2,
        ),
        context_size=VALID.get("context_size", 8),
        context_flag=(
            f"--context-accessions '{VALID['context_accessions']}'"
            if VALID.get("context_accessions")
            else ""
        ),
        occupancy=VALID.get(
            "min_column_occupancy",
            0.70,
        ),
        segment_length=VALID.get(
            "min_segment_length",
            100,
        ),
        model=VALID.get("iqtree_model", "MFP"),
        ufboot=VALID.get("ufboot", 1000),
        alrt=VALID.get("alrt", 1000),
        min_alrt=VALID.get("minimum_alrt", 80),
        min_ufboot=VALID.get("minimum_ufboot", 95),
        formats=VALID.get("formats", "png,pdf"),
        reuse_flag=(
            "--reuse-existing"
            if VALID.get("reuse_existing", False)
            else ""
        ),
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "{params.outdir}" \
            "$(dirname "{log}")" \
            "{params.benchmark_dir}"

        if [ ! -s "{input.events}" ]; then
            echo \
                "ERROR: tableau de concordance absent ou vide : {input.events}" \
                >&2
            exit 1
        fi

        if [ ! -s "{input.alignment}" ]; then
            echo \
                "ERROR: alignement absent ou vide : {input.alignment}" \
                >&2
            exit 1
        fi

        if [ ! -s "{input.script}" ]; then
            echo \
                "ERROR: script de validation absent ou vide : {input.script}" \
                >&2
            exit 1
        fi

        mplconfig=$(mktemp -d \
            "${{TMPDIR:-/tmp}}/matplotlib-vipan-validation.XXXXXX")

        cleanup() {{
            rm -rf "$mplconfig"
        }}

        trap cleanup EXIT

        export MPLCONFIGDIR="$mplconfig"

        python "{input.script}" \
            --events "{input.events}" \
            --alignment "{input.alignment}" \
            --outdir "{params.outdir}" \
            --tiers "{params.tiers}" \
            --topology "{params.topology}" \
            --window {params.window} \
            --step {params.step} \
            --bootstrap-replicates {params.bootstrap} \
            --processes {threads} \
            --seed {params.seed} \
            --rdp-candidates-per-role {params.rdp_candidates} \
            --context-size {params.context_size} \
            --min-column-occupancy {params.occupancy} \
            --min-segment-length {params.segment_length} \
            --formats "{params.formats}" \
            --run-iqtree \
            --iqtree-model "{params.model}" \
            --ufboot {params.ufboot} \
            --alrt {params.alrt} \
            --iqtree-threads {threads} \
            --validation-min-alrt {params.min_alrt} \
            --validation-min-ufboot {params.min_ufboot} \
            {params.context_flag} \
            {params.reuse_flag} \
            > "{log}" 2>&1

        if [ ! -s "{output.manifest}" ]; then
            echo \
                "ERROR: validation_manifest.tsv n'a pas été produit." \
                >&2
            tail -n 100 "{log}" >&2
            exit 1
        fi

        if [ ! -s "{output.tree_summary}" ]; then
            echo \
                "ERROR: segment_tree_summary.tsv n'a pas été produit." \
                >&2
            tail -n 100 "{log}" >&2
            exit 1
        fi

        if [ ! -s "{output.calls}" ]; then
            echo \
                "ERROR: event_validation_calls.tsv n'a pas été produit." \
                >&2
            tail -n 100 "{log}" >&2
            exit 1
        fi
        """

rule classify_validation_calls:
    input:
        calls=f"{RUN_DIR}/validation/event_validation_calls.tsv",
        tree_summary=f"{RUN_DIR}/validation/segment_tree_summary.tsv",
        script="scripts/classify_validation_calls.py",
    output:
        detailed=f"{RUN_DIR}/validation/event_validation_calls_strict.tsv",
        events=f"{RUN_DIR}/validation/event_validation_summary_strict.tsv",
    log:
        f"{RUN_DIR}/logs/classify_validation_calls.log",
    conda:
        "envs/recombination_core.yaml"
    params:
        min_alrt=VALID.get("minimum_alrt", 80),
        min_ufboot=VALID.get("minimum_ufboot", 95),
        local=VALID.get("local_bipartition_maximum", 4),
        intermediate=VALID.get("intermediate_bipartition_maximum", 6),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.detailed}")" "$(dirname "{log}")"
        python "{input.script}" \
            --calls "{input.calls}" \
            --tree-summary "{input.tree_summary}" \
            --output-calls "{output.detailed}" \
            --output-events "{output.events}" \
            --minimum-alrt {params.min_alrt} \
            --minimum-ufboot {params.min_ufboot} \
            --local-bipartition-maximum {params.local} \
            --intermediate-bipartition-maximum {params.intermediate} \
            > "{log}" 2>&1
        """


rule analyze_segment_phylogenetic_signal:
    input:
        tree_summary=f"{RUN_DIR}/validation/segment_tree_summary.tsv",
        strict_calls=f"{RUN_DIR}/validation/event_validation_calls_strict.tsv",
        script="scripts/analyze_segment_phylogenetic_signal.py",
    output:
        segments=f"{RUN_DIR}/validation/phylogenetic_signal/segment_phylogenetic_signal.tsv",
        events=f"{RUN_DIR}/validation/phylogenetic_signal/event_phylogenetic_signal_summary.tsv",
        png=f"{RUN_DIR}/validation/phylogenetic_signal/phylogenetic_signal_diagnostics.png",
        pdf=f"{RUN_DIR}/validation/phylogenetic_signal/phylogenetic_signal_diagnostics.pdf",
    log:
        f"{RUN_DIR}/logs/analyze_segment_phylogenetic_signal.log",
    conda:
        "envs/recombination_validation.yaml"
    params:
        outdir=f"{RUN_DIR}/validation/phylogenetic_signal",
        classes=SIGNAL.get("target_classes", "D_expected_switch_insufficient_branch_support"),
        occupancy=SIGNAL.get("minimum_site_occupancy", 0.70),
        informative=SIGNAL.get("minimum_informative_sites", 20),
        informative_percent=SIGNAL.get("minimum_informative_percent", 1.0),
        comparable=SIGNAL.get("minimum_comparable_sites", 50),
        margin=SIGNAL.get("equivalent_identity_margin", 0.5),
        min_alrt=VALID.get("minimum_alrt", 80),
        min_ufboot=VALID.get("minimum_ufboot", 95),
        formats=SIGNAL.get("formats", "png,pdf"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        export MPLCONFIGDIR="${{TMPDIR:-/tmp}}/matplotlib-vipan-signal"
        python "{input.script}" \
            --tree-summary "{input.tree_summary}" \
            --strict-calls "{input.strict_calls}" \
            --output-dir "{params.outdir}" \
            --target-classes "{params.classes}" \
            --minimum-site-occupancy {params.occupancy} \
            --minimum-informative-sites {params.informative} \
            --minimum-informative-percent {params.informative_percent} \
            --minimum-comparable-sites {params.comparable} \
            --equivalent-identity-margin {params.margin} \
            --minimum-alrt {params.min_alrt} \
            --minimum-ufboot {params.min_ufboot} \
            --formats "{params.formats}" \
            > "{log}" 2>&1
        """


rule write_run_manifest:
    input:
        config="config.yaml",
        raw=config["input"]["fasta"],
        final_alignment=f"{RUN_DIR}/alignment/final_alignment.fasta",
        filter_decisions=f"{RUN_DIR}/qc/distance_filter_decisions.tsv",
        alignment_summary=f"{RUN_DIR}/qc/final_alignment_summary.tsv",
        threeseq=f"{RUN_DIR}/recombination/3seq_output.rec",
        rdp=f"{RUN_DIR}/recombination/rdp_events_robust.tsv",
        concordance=f"{RUN_DIR}/recombination/concordance.summary.tsv",
        strict=f"{RUN_DIR}/validation/event_validation_summary_strict.tsv",
        signal=f"{RUN_DIR}/validation/phylogenetic_signal/event_phylogenetic_signal_summary.tsv",
        executable=THREESEQ["executable"],
        ptable=THREESEQ["ptable"],
        script="scripts/write_run_manifest.py",
    output:
        f"{RUN_DIR}/provenance/run_manifest.json",
    log:
        f"{RUN_DIR}/logs/write_run_manifest.log",
    conda:
        "envs/recombination_core.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"
        python "{input.script}" \
            --output "{output}" \
            --run-id "{RUN_ID}" \
            --file "config={input.config}" \
            --file "raw_fasta={input.raw}" \
            --file "final_alignment={input.final_alignment}" \
            --file "distance_filter_decisions={input.filter_decisions}" \
            --file "alignment_summary={input.alignment_summary}" \
            --file "3seq_output={input.threeseq}" \
            --file "rdp_robust_events={input.rdp}" \
            --file "concordance_summary={input.concordance}" \
            --file "strict_validation_summary={input.strict}" \
            --file "phylogenetic_signal_summary={input.signal}" \
            --file "3seq_executable={input.executable}" \
            --file "3seq_ptable={input.ptable}" \
            --version-command "python=python --version" \
            --version-command "mafft=mafft --version" \
            > "{log}" 2>&1
        """

include: "workflow/phylogenetic_network.smk"
