#!/usr/bin/env python3
"""Compare medoid- and reference-anchored alignment-filtering decisions.

The comparison is accession matched. It reports agreement, discordant
decisions, Jaccard similarity of the retained sets, Cohen's kappa, and a
diagnostic figure. It never changes either source analysis.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-vipan-comparison")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REQUIRED_FIELDS = {
    "accession",
    "status",
    "reason",
    "anchor_accession",
    "distance_to_anchor",
    "comparable_sites_to_anchor",
}


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def as_float(value: object) -> Optional[float]:
    try:
        result = float(clean(value))
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def read_decisions(path: Path) -> Dict[str, Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = REQUIRED_FIELDS - fields
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )
        rows: Dict[str, Dict[str, str]] = {}
        for row in reader:
            accession = clean(row.get("accession"))
            if not accession:
                raise ValueError(f"empty accession in {path}")
            if accession in rows:
                raise ValueError(f"duplicate accession {accession} in {path}")
            rows[accession] = dict(row)
    if not rows:
        raise ValueError(f"no filtering decision in {path}")
    return rows


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def decision_class(medoid_status: str, reference_status: str) -> str:
    if medoid_status == "kept" and reference_status == "kept":
        return "kept_by_both"
    if medoid_status == "excluded" and reference_status == "excluded":
        return "excluded_by_both"
    if medoid_status == "kept" and reference_status == "excluded":
        return "medoid_only"
    if medoid_status == "excluded" and reference_status == "kept":
        return "reference_only"
    return "missing_or_unclassified"


def cohen_kappa(rows: Sequence[Mapping[str, object]]) -> Optional[float]:
    common = [
        row
        for row in rows
        if row["Medoid_status"] in {"kept", "excluded"}
        and row["Reference_status"] in {"kept", "excluded"}
    ]
    if not common:
        return None
    total = len(common)
    observed = sum(
        row["Medoid_status"] == row["Reference_status"] for row in common
    ) / total
    medoid_kept = sum(row["Medoid_status"] == "kept" for row in common) / total
    reference_kept = sum(row["Reference_status"] == "kept" for row in common) / total
    expected = (
        medoid_kept * reference_kept
        + (1.0 - medoid_kept) * (1.0 - reference_kept)
    )
    if math.isclose(expected, 1.0):
        return None
    return (observed - expected) / (1.0 - expected)


def save_figure(
    rows: Sequence[Mapping[str, object]],
    counts: Counter[str],
    output_prefix: Path,
    medoid_threshold: float,
    reference_threshold: float,
    formats: Iterable[str],
) -> None:
    colors = {
        "kept_by_both": "#2ca02c",
        "excluded_by_both": "#7f7f7f",
        "medoid_only": "#1f77b4",
        "reference_only": "#d62728",
        "missing_or_unclassified": "#9467bd",
    }
    labels = {
        "kept_by_both": "Kept by both",
        "excluded_by_both": "Excluded by both",
        "medoid_only": "Medoid only",
        "reference_only": "Reference only",
        "missing_or_unclassified": "Missing/unclassified",
    }
    figure, axes = plt.subplots(1, 2, figsize=(14, 6))

    plotted = 0
    for category in labels:
        selected = [
            row
            for row in rows
            if row["Decision_concordance"] == category
            and isinstance(row["Medoid_distance"], float)
            and isinstance(row["Reference_distance"], float)
        ]
        if not selected:
            continue
        plotted += len(selected)
        axes[0].scatter(
            [row["Medoid_distance"] for row in selected],
            [row["Reference_distance"] for row in selected],
            s=28,
            alpha=0.75,
            label=labels[category],
            color=colors[category],
        )
    axes[0].axvline(medoid_threshold, linestyle="--", linewidth=1.2, color="#333333")
    axes[0].axhline(reference_threshold, linestyle="--", linewidth=1.2, color="#333333")
    axes[0].set_xlabel("p-distance to medoid")
    axes[0].set_ylabel("p-distance to NC_004650.1")
    axes[0].set_title(f"Matched distance profiles (n={plotted})")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8, loc="best")

    order = [
        "kept_by_both",
        "excluded_by_both",
        "medoid_only",
        "reference_only",
        "missing_or_unclassified",
    ]
    present = [category for category in order if counts.get(category, 0)]
    axes[1].barh(
        [labels[category] for category in present][::-1],
        [counts[category] for category in present][::-1],
        color=[colors[category] for category in present][::-1],
    )
    axes[1].set_xlabel("Number of accessions")
    axes[1].set_title("Filtering-decision concordance")
    axes[1].grid(axis="x", alpha=0.2)

    figure.suptitle("Sensitivity to alignment-filtering anchor")
    figure.tight_layout()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        figure.savefig(
            f"{output_prefix}.{extension}",
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
        )
    plt.close(figure)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare medoid- and reference-anchored sequence filtering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--medoid", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-accession", default="NC_004650.1")
    parser.add_argument("--medoid-threshold", type=float, default=0.25)
    parser.add_argument("--reference-threshold", type=float, default=0.25)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--formats", default="png,pdf")
    args = parser.parse_args(argv)
    for name in ("medoid_threshold", "reference_threshold"):
        value = getattr(args, name)
        if not 0.0 <= value <= 1.0:
            parser.error(f"--{name.replace('_', '-')} must be in [0,1]")
    args.formats = tuple(
        item.strip().lower() for item in args.formats.split(",") if item.strip()
    )
    unsupported = set(args.formats) - {"png", "pdf"}
    if unsupported:
        parser.error(f"unsupported format(s): {', '.join(sorted(unsupported))}")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        medoid = read_decisions(args.medoid)
        reference = read_decisions(args.reference)
        all_accessions = sorted(set(medoid) | set(reference))
        rows = []
        for accession in all_accessions:
            medoid_row = medoid.get(accession, {})
            reference_row = reference.get(accession, {})
            medoid_status = clean(medoid_row.get("status")) or "not_present"
            reference_status = clean(reference_row.get("status")) or "not_present"
            category = decision_class(medoid_status, reference_status)
            rows.append(
                {
                    "Accession": accession,
                    "Medoid_status": medoid_status,
                    "Medoid_reason": clean(medoid_row.get("reason")),
                    "Medoid_anchor": clean(medoid_row.get("anchor_accession")),
                    "Medoid_distance": as_float(medoid_row.get("distance_to_anchor")),
                    "Medoid_comparable_sites": clean(
                        medoid_row.get("comparable_sites_to_anchor")
                    ),
                    "Reference_status": reference_status,
                    "Reference_reason": clean(reference_row.get("reason")),
                    "Reference_anchor": clean(reference_row.get("anchor_accession")),
                    "Reference_distance": as_float(
                        reference_row.get("distance_to_anchor")
                    ),
                    "Reference_comparable_sites": clean(
                        reference_row.get("comparable_sites_to_anchor")
                    ),
                    "Decision_concordance": category,
                }
            )

        fields = [
            "Accession",
            "Medoid_status",
            "Medoid_reason",
            "Medoid_anchor",
            "Medoid_distance",
            "Medoid_comparable_sites",
            "Reference_status",
            "Reference_reason",
            "Reference_anchor",
            "Reference_distance",
            "Reference_comparable_sites",
            "Decision_concordance",
        ]
        serialized = []
        for row in rows:
            current = dict(row)
            for field in ("Medoid_distance", "Reference_distance"):
                if isinstance(current[field], float):
                    current[field] = f"{current[field]:.8f}"
                elif current[field] is None:
                    current[field] = ""
            serialized.append(current)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(
            args.output_dir / "sequence_filtering_comparison.tsv",
            serialized,
            fields,
        )
        discordant = [
            row
            for row in serialized
            if row["Decision_concordance"]
            in {"medoid_only", "reference_only", "missing_or_unclassified"}
        ]
        write_tsv(
            args.output_dir / "discordant_accessions.tsv",
            discordant,
            fields,
        )

        counts = Counter(clean(row["Decision_concordance"]) for row in rows)
        medoid_kept = {key for key, row in medoid.items() if clean(row.get("status")) == "kept"}
        reference_kept = {
            key for key, row in reference.items() if clean(row.get("status")) == "kept"
        }
        kept_union = medoid_kept | reference_kept
        kept_intersection = medoid_kept & reference_kept
        common = set(medoid) & set(reference)
        agreement = (
            sum(
                clean(medoid[key].get("status"))
                == clean(reference[key].get("status"))
                for key in common
            )
            / len(common)
            if common
            else float("nan")
        )
        kappa = cohen_kappa(rows)
        summary_rows = [
            {"Metric": "reference_accession", "Value": args.reference_accession},
            {"Metric": "accessions_union", "Value": len(all_accessions)},
            {"Metric": "accessions_common", "Value": len(common)},
            {"Metric": "medoid_kept", "Value": len(medoid_kept)},
            {"Metric": "reference_kept", "Value": len(reference_kept)},
            {"Metric": "kept_by_both", "Value": counts["kept_by_both"]},
            {"Metric": "excluded_by_both", "Value": counts["excluded_by_both"]},
            {"Metric": "medoid_only", "Value": counts["medoid_only"]},
            {"Metric": "reference_only", "Value": counts["reference_only"]},
            {
                "Metric": "decision_agreement_percent",
                "Value": "" if math.isnan(agreement) else f"{100.0 * agreement:.6f}",
            },
            {
                "Metric": "kept_set_jaccard",
                "Value": (
                    f"{len(kept_intersection) / len(kept_union):.8f}"
                    if kept_union
                    else ""
                ),
            },
            {
                "Metric": "cohen_kappa",
                "Value": "" if kappa is None else f"{kappa:.8f}",
            },
        ]
        write_tsv(
            args.output_dir / "filtering_comparison_summary.tsv",
            summary_rows,
            ["Metric", "Value"],
        )
        save_figure(
            rows,
            counts,
            args.output_dir / "filtering_strategy_comparison",
            args.medoid_threshold,
            args.reference_threshold,
            args.formats,
        )

        print(f"Accessions in union       : {len(all_accessions)}")
        print(f"Kept by medoid            : {len(medoid_kept)}")
        print(f"Kept by reference         : {len(reference_kept)}")
        print(f"Kept by both              : {counts['kept_by_both']}")
        print(f"Medoid only               : {counts['medoid_only']}")
        print(f"Reference only            : {counts['reference_only']}")
        print(f"Decision agreement        : {100.0 * agreement:.3f}%")
        print(
            "Kept-set Jaccard          : "
            + (
                f"{len(kept_intersection) / len(kept_union):.6f}"
                if kept_union
                else "not defined"
            )
        )
        print(
            "Cohen's kappa             : "
            + (f"{kappa:.6f}" if kappa is not None else "not defined")
        )
        print(f"Output directory          : {args.output_dir}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
