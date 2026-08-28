"""Compare completed strict validation results for both filtering anchors."""


MEDOID_RUN = "splcv_refiltered_medoid_deterministic_v2"
REFERENCE_RUN = "splcv_refiltered_reference_NC004650_deterministic_v2"
OUTPUT_DIR = "results/sensitivity_medoid_vs_NC004650_deterministic_v2"


rule all:
    input:
        f"{OUTPUT_DIR}/strict_validation_comparison.tsv",
        f"{OUTPUT_DIR}/strict_validation_comparison_summary.tsv",
        f"{OUTPUT_DIR}/strict_validation_class_comparison.png",
        f"{OUTPUT_DIR}/strict_validation_class_comparison.pdf",


rule compare_validation_strategies:
    input:
        medoid=(
            f"results/{MEDOID_RUN}/validation/"
            "event_validation_summary_strict.tsv"
        ),
        reference=(
            f"results/{REFERENCE_RUN}/validation/"
            "event_validation_summary_strict.tsv"
        ),
        script="scripts/compare_validation_strategies.py",
    output:
        detailed=f"{OUTPUT_DIR}/strict_validation_comparison.tsv",
        summary=f"{OUTPUT_DIR}/strict_validation_comparison_summary.tsv",
        png=f"{OUTPUT_DIR}/strict_validation_class_comparison.png",
        pdf=f"{OUTPUT_DIR}/strict_validation_class_comparison.pdf",
    log:
        f"{OUTPUT_DIR}/compare_validation_strategies.log",
    conda:
        "../envs/recombination_core.yaml"
    params:
        outdir=OUTPUT_DIR,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        export MPLCONFIGDIR="${{TMPDIR:-/tmp}}/matplotlib-vipan-validation-sensitivity"

        python "{input.script}" \
            --medoid "{input.medoid}" \
            --reference "{input.reference}" \
            --output-dir "{params.outdir}" \
            > "{log}" 2>&1

        test -s "{output.detailed}"
        test -s "{output.summary}"
        test -s "{output.png}"
        test -s "{output.pdf}"
        """
