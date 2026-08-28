#!/usr/bin/env python3
"""Filter distant or poor-quality sequences from a nucleotide alignment.

Distances are uncorrected p-distances calculated only at A/C/G/T positions.
Every decision is written to a TSV audit table.  The filtered FASTA is emitted
without alignment gaps so it can be realigned from scratch downstream.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-vipan-filter")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fasta_utils import FastaRecord, read_fasta, write_fasta


VALID_BASES = frozenset("ACGT")
GAPS = frozenset("-.")


@dataclass(frozen=True)
class Distance:
    p_distance: Optional[float]
    comparable_sites: int


def first_token(value: str) -> str:
    return str(value).strip().split()[0] if str(value).strip() else ""


def read_id_file(path: Optional[Path]) -> Set[str]:
    if path is None:
        return set()
    identifiers = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            value = raw_line.strip()
            if value and not value.startswith("#"):
                identifiers.add(first_token(value))
    return identifiers


def distance(sequence_a: str, sequence_b: str) -> Distance:
    comparable = 0
    differences = 0
    for base_a, base_b in zip(sequence_a, sequence_b):
        if base_a in VALID_BASES and base_b in VALID_BASES:
            comparable += 1
            differences += base_a != base_b
    return Distance(differences / comparable if comparable else None, comparable)


def pair_key(first: str, second: str) -> Tuple[str, str]:
    return (first, second) if first < second else (second, first)


def distance_lookup(
    first: str,
    second: str,
    distances: Mapping[Tuple[str, str], Distance],
) -> Distance:
    if first == second:
        return Distance(0.0, 0)
    return distances[pair_key(first, second)]


def select_medoid(
    candidates: Sequence[str],
    distances: Mapping[Tuple[str, str], Distance],
    minimum_comparable: int,
) -> Tuple[str, Dict[str, float]]:
    medians: Dict[str, float] = {}
    for accession in candidates:
        values = []
        for other in candidates:
            if other == accession:
                continue
            current = distance_lookup(accession, other, distances)
            if current.p_distance is not None and current.comparable_sites >= minimum_comparable:
                values.append(current.p_distance)
        if values:
            medians[accession] = statistics.median(values)
    if not medians:
        raise ValueError("no medoid can be calculated with the current overlap threshold")
    return min(medians, key=lambda item: (medians[item], item)), medians


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove alignment outliers with a fully audited distance rule.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--kept-fasta", required=True, type=Path)
    parser.add_argument("--kept-aligned-fasta", required=True, type=Path)
    parser.add_argument("--excluded-fasta", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--pairwise-distances", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--plot-prefix", required=True, type=Path)
    parser.add_argument("--mode", choices=("medoid", "reference", "none"), default="medoid")
    parser.add_argument("--reference-accession")
    parser.add_argument("--max-p-distance", type=float, default=0.25)
    parser.add_argument("--min-comparable-sites", type=int, default=1000)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.60)
    parser.add_argument("--min-valid-bases", type=int, default=2400)
    parser.add_argument("--max-gap-fraction", type=float, default=0.20)
    parser.add_argument("--max-ambiguity-fraction", type=float, default=0.05)
    parser.add_argument("--manual-keep", type=Path)
    parser.add_argument("--manual-exclude", type=Path)
    args = parser.parse_args(argv)
    if not 0 <= args.max_p_distance <= 1:
        parser.error("--max-p-distance must be in [0,1]")
    for name in ("min_overlap_fraction", "max_gap_fraction", "max_ambiguity_fraction"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0,1]")
    if args.min_comparable_sites < 1 or args.min_valid_bases < 1:
        parser.error("site-count thresholds must be positive")
    if args.mode == "reference" and not args.reference_accession:
        parser.error("--reference-accession is required when --mode reference")
    return args


def write_plot(
    rows: Sequence[Mapping[str, object]],
    maximum_distance: float,
    prefix: Path,
) -> None:
    usable = [row for row in rows if isinstance(row["distance_to_anchor"], float)]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    kept = [row for row in usable if row["status"] == "kept"]
    excluded = [row for row in usable if row["status"] != "kept"]

    axes[0].hist(
        [row["distance_to_anchor"] for row in usable],
        bins=min(40, max(8, int(math.sqrt(max(1, len(usable)))) * 2)),
        color="#4c78a8",
        edgecolor="white",
    )
    axes[0].axvline(maximum_distance, color="#d62728", linestyle="--", label="configured threshold")
    axes[0].set_xlabel("p-distance to anchor")
    axes[0].set_ylabel("Number of sequences")
    axes[0].legend()

    for subset, color, label in (
        (kept, "#2ca02c", "kept"),
        (excluded, "#d62728", "excluded"),
    ):
        axes[1].scatter(
            [row["distance_to_anchor"] for row in subset],
            [row["gap_fraction"] for row in subset],
            s=25,
            alpha=0.75,
            color=color,
            label=label,
        )
    axes[1].axvline(maximum_distance, color="#333333", linestyle="--", linewidth=1)
    axes[1].set_xlabel("p-distance to anchor")
    axes[1].set_ylabel("Gap fraction in initial alignment")
    axes[1].legend()
    axes[1].grid(alpha=0.2)
    figure.tight_layout()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    figure.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(figure)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        records = read_fasta(args.alignment)
        if len(records) < 3:
            raise ValueError("at least three aligned sequences are required")
        lengths = {len(record.sequence) for record in records}
        if len(lengths) != 1:
            raise ValueError(f"input FASTA is not aligned; lengths={sorted(lengths)[:10]}")
        alignment_length = lengths.pop()
        minimum_comparable = max(
            args.min_comparable_sites,
            math.ceil(args.min_overlap_fraction * alignment_length),
        )

        sequences: Dict[str, str] = {}
        descriptions: Dict[str, str] = {}
        base_metrics: Dict[str, Dict[str, object]] = {}
        for record in records:
            accession = first_token(record.id)
            if accession in sequences:
                raise ValueError(f"duplicate alignment identifier: {accession}")
            sequence = record.sequence.upper().replace("U", "T")
            sequences[accession] = sequence
            descriptions[accession] = record.description
            valid = sum(base in VALID_BASES for base in sequence)
            gaps = sum(base in GAPS for base in sequence)
            ambiguous = alignment_length - valid - gaps
            base_metrics[accession] = {
                "valid_bases": valid,
                "gap_fraction": gaps / alignment_length,
                "ambiguity_fraction": ambiguous / alignment_length,
            }

        manual_keep = read_id_file(args.manual_keep)
        manual_exclude = read_id_file(args.manual_exclude)
        unknown_manual = (manual_keep | manual_exclude) - set(sequences)
        if unknown_manual:
            raise ValueError(
                "manual include/exclude identifiers absent from alignment: "
                + ", ".join(sorted(unknown_manual)[:20])
            )
        if manual_keep & manual_exclude:
            raise ValueError("some accessions occur in both manual keep and exclude lists")

        base_pass = []
        preliminary_reasons: Dict[str, str] = {}
        for accession, metrics in base_metrics.items():
            if metrics["valid_bases"] < args.min_valid_bases:
                preliminary_reasons[accession] = "too_few_valid_bases"
            elif metrics["gap_fraction"] > args.max_gap_fraction:
                preliminary_reasons[accession] = "excessive_gap_fraction"
            elif metrics["ambiguity_fraction"] > args.max_ambiguity_fraction:
                preliminary_reasons[accession] = "excessive_ambiguity_fraction"
            elif accession in manual_exclude:
                preliminary_reasons[accession] = "manual_exclusion"
            else:
                base_pass.append(accession)
        if len(base_pass) < 3:
            raise ValueError("fewer than three sequences pass basic alignment QC")

        pairwise: Dict[Tuple[str, str], Distance] = {}
        pair_rows = []
        accessions = sorted(sequences)
        for index, first in enumerate(accessions):
            for second in accessions[index + 1 :]:
                current = distance(sequences[first], sequences[second])
                pairwise[(first, second)] = current
                pair_rows.append(
                    {
                        "accession_1": first,
                        "accession_2": second,
                        "comparable_sites": current.comparable_sites,
                        "p_distance": "" if current.p_distance is None else f"{current.p_distance:.8f}",
                        "identity_percent": "" if current.p_distance is None else f"{100 * (1-current.p_distance):.6f}",
                    }
                )

        medoid_medians: Dict[str, float] = {}
        if args.mode == "reference":
            anchor = first_token(args.reference_accession)
            if anchor not in base_pass:
                raise ValueError(f"reference accession {anchor} did not pass basic QC")
        elif args.mode == "medoid":
            anchor, medoid_medians = select_medoid(base_pass, pairwise, minimum_comparable)
        else:
            anchor = ""

        rows = []
        kept_ids = []
        for accession in accessions:
            metrics = base_metrics[accession]
            status = "kept"
            reason = "passed_all_filters"
            anchor_distance: Optional[float] = 0.0 if accession == anchor else None
            anchor_comparable = alignment_length if accession == anchor else 0

            if accession in preliminary_reasons:
                status = "excluded"
                reason = preliminary_reasons[accession]
            elif args.mode != "none" and accession != anchor:
                current = distance_lookup(accession, anchor, pairwise)
                anchor_distance = current.p_distance
                anchor_comparable = current.comparable_sites
                if current.comparable_sites < minimum_comparable or current.p_distance is None:
                    status = "excluded"
                    reason = "insufficient_overlap_with_anchor"
                elif current.p_distance > args.max_p_distance and accession not in manual_keep:
                    status = "excluded"
                    reason = "too_distant_from_anchor"
                elif accession in manual_keep:
                    reason = "manual_keep_override_distance"
            elif accession in manual_keep:
                reason = "manual_keep"

            if status == "kept":
                kept_ids.append(accession)
            rows.append(
                {
                    "accession": accession,
                    "alignment_length": alignment_length,
                    "valid_bases": metrics["valid_bases"],
                    "gap_fraction": metrics["gap_fraction"],
                    "ambiguity_fraction": metrics["ambiguity_fraction"],
                    "filter_mode": args.mode,
                    "anchor_accession": anchor,
                    "distance_to_anchor": anchor_distance,
                    "comparable_sites_to_anchor": anchor_comparable,
                    "median_distance_to_dataset": medoid_medians.get(accession, ""),
                    "status": status,
                    "reason": reason,
                }
            )

        if len(kept_ids) < 3:
            raise ValueError(f"only {len(kept_ids)} sequences remain after distance filtering")

        kept_set = set(kept_ids)
        kept_aligned = []
        kept_unaligned = []
        excluded = []
        for accession in accessions:
            aligned = sequences[accession]
            aligned_record = FastaRecord(accession, accession, aligned)
            if accession in kept_set:
                kept_aligned.append(aligned_record)
                unaligned = aligned.replace("-", "").replace(".", "")
                kept_unaligned.append(FastaRecord(accession, accession, unaligned))
            else:
                excluded.append(aligned_record)

        output_paths = (
            args.kept_fasta,
            args.kept_aligned_fasta,
            args.excluded_fasta,
            args.metrics,
            args.pairwise_distances,
            args.summary_json,
        )
        for path in output_paths:
            path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(kept_unaligned, args.kept_fasta)
        write_fasta(kept_aligned, args.kept_aligned_fasta)
        write_fasta(excluded, args.excluded_fasta)

        metrics_fields = list(rows[0])
        with args.metrics.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=metrics_fields, delimiter="\t")
            writer.writeheader()
            for row in rows:
                serialized = dict(row)
                for field in ("gap_fraction", "ambiguity_fraction", "distance_to_anchor", "median_distance_to_dataset"):
                    if isinstance(serialized[field], float):
                        serialized[field] = f"{serialized[field]:.8f}"
                writer.writerow(serialized)
        with args.pairwise_distances.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(pair_rows)

        reason_counts: Dict[str, int] = {}
        for row in rows:
            reason_counts[str(row["reason"])] = reason_counts.get(str(row["reason"]), 0) + 1
        summary = {
            "alignment": str(args.alignment),
            "alignment_length": alignment_length,
            "input_sequences": len(records),
            "kept_sequences": len(kept_ids),
            "excluded_sequences": len(records) - len(kept_ids),
            "filter_mode": args.mode,
            "anchor_accession": anchor,
            "max_p_distance": args.max_p_distance,
            "minimum_comparable_sites_effective": minimum_comparable,
            "min_valid_bases": args.min_valid_bases,
            "max_gap_fraction": args.max_gap_fraction,
            "max_ambiguity_fraction": args.max_ambiguity_fraction,
            "reason_counts": reason_counts,
        }
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        write_plot(rows, args.max_p_distance, args.plot_prefix)

        print(f"Filter mode       : {args.mode}")
        print(f"Anchor accession  : {anchor or 'not applicable'}")
        print(f"Input sequences   : {len(records)}")
        print(f"Kept sequences    : {len(kept_ids)}")
        print(f"Excluded sequences: {len(records) - len(kept_ids)}")
        print(f"Audit table       : {args.metrics}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
