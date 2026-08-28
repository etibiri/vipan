#!/usr/bin/env python3
"""Compare strict recombination-validation results between two filtering arms."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-vipan-validation-comparison")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CLASS_FIELD = "Best_strict_validation_class"
ACCESSION_FIELD = "Recombinant_accession"
CLASS_ORDER = [
    "A_local_supported_parent_switch",
    "B_intermediate_supported_parent_switch",
    "C_broad_lineage_supported_switch",
    "D_expected_switch_insufficient_branch_support",
    "E_parent_switch_without_expected_roles",
    "F_no_parent_switch_or_unavailable",
]
SUPPORTED_CLASSES = set(CLASS_ORDER[:3])


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def read_summary(path: Path) -> Dict[str, Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {ACCESSION_FIELD, CLASS_FIELD}
        if not required <= fields:
            raise ValueError(
                f"{path} is missing: {', '.join(sorted(required - fields))}"
            )
        rows: Dict[str, Dict[str, str]] = {}
        for row in reader:
            accession = clean(row.get(ACCESSION_FIELD))
            if not accession:
                raise ValueError(f"empty recombinant accession in {path}")
            if accession in rows:
                raise ValueError(f"duplicate recombinant {accession} in {path}")
            rows[accession] = dict(row)
    return rows


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_counts(
    medoid_counts: Counter[str],
    reference_counts: Counter[str],
    prefix: Path,
) -> None:
    labels = [item.split("_", 1)[0] for item in CLASS_ORDER]
    positions = list(range(len(CLASS_ORDER)))
    width = 0.38
    figure, axis = plt.subplots(figsize=(10, 6))
    axis.bar(
        [value - width / 2 for value in positions],
        [medoid_counts[label] for label in CLASS_ORDER],
        width=width,
        label="Medoid",
        color="#4c78a8",
    )
    axis.bar(
        [value + width / 2 for value in positions],
        [reference_counts[label] for label in CLASS_ORDER],
        width=width,
        label="NC_004650.1 reference",
        color="#f58518",
    )
    axis.set_xticks(positions, labels)
    axis.set_xlabel("Strict validation class")
    axis.set_ylabel("Number of recombinants")
    axis.set_title("Sensitivity of strict recombination validation to filtering anchor")
    axis.grid(axis="y", alpha=0.2)
    axis.legend()
    figure.tight_layout()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    figure.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(figure)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare strict validation classes across filtering strategies."
    )
    parser.add_argument("--medoid", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        medoid = read_summary(args.medoid)
        reference = read_summary(args.reference)
        accessions = sorted(set(medoid) | set(reference))
        rows = []
        for accession in accessions:
            medoid_row = medoid.get(accession, {})
            reference_row = reference.get(accession, {})
            medoid_class = clean(medoid_row.get(CLASS_FIELD)) or "not_evaluated"
            reference_class = clean(reference_row.get(CLASS_FIELD)) or "not_evaluated"
            if accession not in medoid:
                concordance = "reference_only"
            elif accession not in reference:
                concordance = "medoid_only"
            elif medoid_class == reference_class:
                concordance = "same_class"
            else:
                concordance = "different_class"
            rows.append(
                {
                    "Recombinant_accession": accession,
                    "Medoid_class": medoid_class,
                    "Medoid_best_scheme": clean(
                        medoid_row.get("Best_breakpoint_scheme")
                    ),
                    "Reference_class": reference_class,
                    "Reference_best_scheme": clean(
                        reference_row.get("Best_breakpoint_scheme")
                    ),
                    "Validation_concordance": concordance,
                }
            )

        medoid_counts = Counter(
            clean(row.get(CLASS_FIELD)) for row in medoid.values()
        )
        reference_counts = Counter(
            clean(row.get(CLASS_FIELD)) for row in reference.values()
        )
        medoid_supported = {
            accession
            for accession, row in medoid.items()
            if clean(row.get(CLASS_FIELD)) in SUPPORTED_CLASSES
        }
        reference_supported = {
            accession
            for accession, row in reference.items()
            if clean(row.get(CLASS_FIELD)) in SUPPORTED_CLASSES
        }
        supported_union = medoid_supported | reference_supported
        supported_intersection = medoid_supported & reference_supported
        concordance_counts = Counter(
            clean(row["Validation_concordance"]) for row in rows
        )

        args.output_dir.mkdir(parents=True, exist_ok=True)
        detailed_fields = list(rows[0]) if rows else [
            "Recombinant_accession",
            "Medoid_class",
            "Medoid_best_scheme",
            "Reference_class",
            "Reference_best_scheme",
            "Validation_concordance",
        ]
        write_tsv(
            args.output_dir / "strict_validation_comparison.tsv",
            rows,
            detailed_fields,
        )
        summary_rows = [
            {"Metric": "medoid_recombinants_evaluated", "Value": len(medoid)},
            {"Metric": "reference_recombinants_evaluated", "Value": len(reference)},
            {"Metric": "same_class", "Value": concordance_counts["same_class"]},
            {"Metric": "different_class", "Value": concordance_counts["different_class"]},
            {"Metric": "medoid_only", "Value": concordance_counts["medoid_only"]},
            {"Metric": "reference_only", "Value": concordance_counts["reference_only"]},
            {"Metric": "medoid_supported_A_to_C", "Value": len(medoid_supported)},
            {"Metric": "reference_supported_A_to_C", "Value": len(reference_supported)},
            {"Metric": "supported_by_both", "Value": len(supported_intersection)},
            {
                "Metric": "supported_set_jaccard",
                "Value": (
                    f"{len(supported_intersection) / len(supported_union):.8f}"
                    if supported_union
                    else ""
                ),
            },
        ]
        for label in CLASS_ORDER:
            summary_rows.append(
                {
                    "Metric": f"medoid_class_{label.split('_', 1)[0]}",
                    "Value": medoid_counts[label],
                }
            )
            summary_rows.append(
                {
                    "Metric": f"reference_class_{label.split('_', 1)[0]}",
                    "Value": reference_counts[label],
                }
            )
        write_tsv(
            args.output_dir / "strict_validation_comparison_summary.tsv",
            summary_rows,
            ["Metric", "Value"],
        )
        plot_counts(
            medoid_counts,
            reference_counts,
            args.output_dir / "strict_validation_class_comparison",
        )

        print(f"Medoid recombinants evaluated    : {len(medoid)}")
        print(f"Reference recombinants evaluated : {len(reference)}")
        print(f"Same strict class                : {concordance_counts['same_class']}")
        print(f"Different strict class           : {concordance_counts['different_class']}")
        print(f"Medoid supported A-C             : {len(medoid_supported)}")
        print(f"Reference supported A-C          : {len(reference_supported)}")
        print(f"Supported by both                : {len(supported_intersection)}")
        print(f"Output directory                 : {args.output_dir}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
