#!/usr/bin/env python3
"""Diagnose why recombinant segment phylogenies have insufficient support.

The program joins ``event_validation_calls_strict.tsv`` with
``segment_tree_summary.tsv`` and computes, for each targeted segment:

* eligible, constant, variable, singleton-variable and
  parsimony-informative alignment sites;
* mean pairwise nucleotide diversity;
* recombinant-to-parent identity and number of differences;
* the closest competing sequence and the rank of the proposed parent;
* an operational diagnosis separating low phylogenetic signal from parental
  ambiguity or unresolved/conflicting topology.

It does not alter the original validation tables or IQ-TREE results.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(os.environ.get("TMPDIR", "/tmp"), "matplotlib-vipan"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


VALID_BASES = frozenset("ACGT")
DEFAULT_TARGET_CLASS = "D_expected_switch_insufficient_branch_support"


def clean(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def to_float(value: object) -> Optional[float]:
    text = clean(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_int(value: object) -> Optional[int]:
    number = to_float(value)
    return None if number is None else int(number)


def to_bool(value: object) -> bool:
    return clean(value).lower() in {"true", "1", "yes", "y"}


def first_token(value: object) -> str:
    text = clean(value)
    return text.split()[0].upper() if text else ""


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    identifier = ""
    parts: List[str] = []

    def save() -> None:
        if not identifier:
            return
        sequence = "".join(parts).upper()
        if identifier in sequences and sequences[identifier] != sequence:
            raise ValueError(f"Duplicate FASTA identifier with different sequence: {identifier}")
        sequences[identifier] = sequence

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                identifier = first_token(line[1:])
                parts = []
            else:
                parts.append(re.sub(r"\s+", "", line))
        save()

    if not sequences:
        raise ValueError(f"No FASTA sequence in {path}")
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"Unaligned FASTA lengths in {path}: {sorted(lengths)}")
    return sequences


def pairwise_identity(sequence_a: str, sequence_b: str) -> Tuple[Optional[float], int, int]:
    comparable = 0
    differences = 0
    for base_a, base_b in zip(sequence_a, sequence_b):
        if base_a in VALID_BASES and base_b in VALID_BASES:
            comparable += 1
            differences += base_a != base_b
    if comparable == 0:
        return None, 0, 0
    identity = 100.0 * (comparable - differences) / comparable
    return identity, comparable, differences


def mean_pairwise_diversity(sequences: Mapping[str, str]) -> Optional[float]:
    distances: List[float] = []
    for first, second in itertools.combinations(sequences, 2):
        identity, comparable, _ = pairwise_identity(sequences[first], sequences[second])
        if identity is not None and comparable:
            distances.append(1.0 - identity / 100.0)
    return float(np.mean(distances)) if distances else None


def alignment_site_metrics(
    sequences: Mapping[str, str], minimum_occupancy: float
) -> Dict[str, object]:
    names = list(sequences)
    number_sequences = len(names)
    alignment_length = len(next(iter(sequences.values())))
    minimum_valid = max(2, math.ceil(minimum_occupancy * number_sequences))

    eligible = 0
    constant = 0
    variable = 0
    singleton_variable = 0
    informative = 0
    total_entropy = 0.0
    occupancy_values: List[float] = []

    for position in range(alignment_length):
        bases = [
            sequences[name][position]
            for name in names
            if sequences[name][position] in VALID_BASES
        ]
        occupancy_values.append(len(bases) / number_sequences)
        if len(bases) < minimum_valid:
            continue
        eligible += 1
        counts = Counter(bases)
        if len(counts) == 1:
            constant += 1
            continue
        variable += 1
        repeated_states = sum(count >= 2 for count in counts.values())
        if repeated_states >= 2:
            informative += 1
        else:
            singleton_variable += 1
        frequencies = [count / len(bases) for count in counts.values()]
        total_entropy += -sum(frequency * math.log2(frequency) for frequency in frequencies)

    return {
        "Number_sequences": number_sequences,
        "Alignment_length": alignment_length,
        "Minimum_valid_sequences_per_site": minimum_valid,
        "Eligible_sites": eligible,
        "Excluded_low_occupancy_sites": alignment_length - eligible,
        "Mean_site_occupancy_percent": (
            100.0 * float(np.mean(occupancy_values)) if occupancy_values else math.nan
        ),
        "Constant_sites": constant,
        "Variable_sites": variable,
        "Singleton_variable_sites": singleton_variable,
        "Parsimony_informative_sites": informative,
        "Variable_percent_of_eligible": 100.0 * variable / eligible if eligible else math.nan,
        "Informative_percent_of_eligible": (
            100.0 * informative / eligible if eligible else math.nan
        ),
        "Informative_sites_per_100_nt": (
            100.0 * informative / alignment_length if alignment_length else math.nan
        ),
        "Mean_entropy_variable_sites": total_entropy / variable if variable else 0.0,
        "Unique_aligned_haplotypes": len(set(sequences.values())),
        "Mean_pairwise_pi": mean_pairwise_diversity(sequences),
    }


def parent_competition_metrics(
    sequences: Mapping[str, str],
    recombinant: str,
    parent: str,
    equivalent_margin: float,
) -> Dict[str, object]:
    result: Dict[str, object] = {
        "Recombinant_present": recombinant in sequences,
        "Nearest_parent_present": parent in sequences,
        "Recombinant_parent_comparable_sites": "",
        "Recombinant_parent_differences": "",
        "Recombinant_parent_identity_percent": "",
        "Closest_sequence_overall": "",
        "Closest_sequence_identity_percent": "",
        "Best_competing_sequence": "",
        "Best_competing_identity_percent": "",
        "Parent_identity_margin_percentage_points": "",
        "Parent_identity_rank": "",
        "Near_equivalent_sequences_count": "",
        "Near_equivalent_sequences": "",
    }
    if recombinant not in sequences or parent not in sequences:
        return result

    parent_identity, comparable, differences = pairwise_identity(
        sequences[recombinant], sequences[parent]
    )
    result["Recombinant_parent_comparable_sites"] = comparable
    result["Recombinant_parent_differences"] = differences
    result["Recombinant_parent_identity_percent"] = parent_identity

    ranked: List[Tuple[str, float]] = []
    for accession, sequence in sequences.items():
        if accession == recombinant:
            continue
        identity, valid, _ = pairwise_identity(sequences[recombinant], sequence)
        if identity is not None and valid:
            ranked.append((accession, identity))
    ranked.sort(key=lambda item: (-item[1], item[0]))
    if not ranked or parent_identity is None:
        return result

    result["Closest_sequence_overall"] = ranked[0][0]
    result["Closest_sequence_identity_percent"] = ranked[0][1]
    competitors = [item for item in ranked if item[0] != parent]
    if competitors:
        result["Best_competing_sequence"] = competitors[0][0]
        result["Best_competing_identity_percent"] = competitors[0][1]
        result["Parent_identity_margin_percentage_points"] = (
            parent_identity - competitors[0][1]
        )
    result["Parent_identity_rank"] = next(
        index for index, item in enumerate(ranked, start=1) if item[0] == parent
    )
    near_equivalent = [
        accession
        for accession, identity in ranked
        if identity >= parent_identity - equivalent_margin
    ]
    result["Near_equivalent_sequences_count"] = len(near_equivalent)
    result["Near_equivalent_sequences"] = ";".join(near_equivalent)
    return result


def support_passes(
    alrt: Optional[float],
    ufboot: Optional[float],
    minimum_alrt: float,
    minimum_ufboot: float,
) -> bool:
    return (
        alrt is not None
        and ufboot is not None
        and alrt >= minimum_alrt
        and ufboot >= minimum_ufboot
    )


def diagnose_segment(
    metrics: Mapping[str, object],
    minimum_informative_sites: int,
    minimum_informative_percent: float,
    minimum_comparable_sites: int,
    competition_margin: float,
) -> str:
    if to_bool(metrics.get("Branch_support_pass")):
        return "Supported_segment"
    alignment_length = to_int(metrics.get("Alignment_length")) or 0
    eligible = to_int(metrics.get("Eligible_sites")) or 0
    informative = to_int(metrics.get("Parsimony_informative_sites")) or 0
    informative_percent = to_float(metrics.get("Informative_percent_of_eligible")) or 0.0
    comparable = to_int(metrics.get("Recombinant_parent_comparable_sites")) or 0
    parent_rank = to_int(metrics.get("Parent_identity_rank"))
    margin = to_float(metrics.get("Parent_identity_margin_percentage_points"))
    equivalent_count = to_int(metrics.get("Near_equivalent_sequences_count")) or 0

    if alignment_length == 0 or eligible == 0:
        return "Missing_or_unusable_alignment"
    if comparable < minimum_comparable_sites:
        return "Insufficient_recombinant_parent_overlap"
    if informative < minimum_informative_sites or informative_percent < minimum_informative_percent:
        return "Low_phylogenetic_signal"
    if parent_rank is not None and parent_rank > 1 and margin is not None and margin < -competition_margin:
        return "Closer_competing_sequence"
    if equivalent_count >= 2 and margin is not None and abs(margin) <= competition_margin:
        return "Parent_indistinguishable_within_lineage"
    return "Adequate_signal_but_unresolved_or_conflicting_branch"


def key(row: Mapping[str, object], include_segment: bool = False) -> Tuple[str, ...]:
    values = (
        clean(row.get("Recombinant_accession")),
        clean(row.get("3SEQ_Event_ID")),
        clean(row.get("RDP_Event_ID")),
        clean(row.get("Breakpoint_scheme")),
    )
    return (*values, clean(row.get("Segment"))) if include_segment else values


SEGMENT_FIELDS = [
    "Evidence_tier",
    "Recombinant_accession",
    "3SEQ_Event_ID",
    "RDP_Event_ID",
    "Breakpoint_scheme",
    "Segment",
    "Strict_validation_class",
    "Nearest_selected_parent",
    "Nearest_parent_roles",
    "SH_aLRT",
    "UFBoot",
    "Branch_support_pass",
    "Supported_bipartition_size",
    "Tree_number_taxa",
    "Segment_fasta",
    "Number_sequences",
    "Alignment_length",
    "Minimum_valid_sequences_per_site",
    "Eligible_sites",
    "Excluded_low_occupancy_sites",
    "Mean_site_occupancy_percent",
    "Constant_sites",
    "Variable_sites",
    "Singleton_variable_sites",
    "Parsimony_informative_sites",
    "Variable_percent_of_eligible",
    "Informative_percent_of_eligible",
    "Informative_sites_per_100_nt",
    "Mean_entropy_variable_sites",
    "Unique_aligned_haplotypes",
    "Mean_pairwise_pi",
    "Recombinant_present",
    "Nearest_parent_present",
    "Recombinant_parent_comparable_sites",
    "Recombinant_parent_differences",
    "Recombinant_parent_identity_percent",
    "Closest_sequence_overall",
    "Closest_sequence_identity_percent",
    "Best_competing_sequence",
    "Best_competing_identity_percent",
    "Parent_identity_margin_percentage_points",
    "Parent_identity_rank",
    "Near_equivalent_sequences_count",
    "Near_equivalent_sequences",
    "Signal_diagnosis",
]


CALL_FIELDS = [
    "Evidence_tier",
    "Recombinant_accession",
    "3SEQ_Event_ID",
    "RDP_Event_ID",
    "Breakpoint_scheme",
    "Strict_validation_class",
    "Role_switch_evidence",
    "Segment_A_signal_diagnosis",
    "Segment_A_informative_sites",
    "Segment_A_informative_percent",
    "Segment_A_parent_identity_percent",
    "Segment_A_parent_rank",
    "Segment_B_signal_diagnosis",
    "Segment_B_informative_sites",
    "Segment_B_informative_percent",
    "Segment_B_parent_identity_percent",
    "Segment_B_parent_rank",
    "Combined_signal_diagnosis",
]


def combined_diagnosis(segment_a: Mapping[str, object], segment_b: Mapping[str, object]) -> str:
    diagnoses = {
        clean(segment_a.get("Signal_diagnosis")),
        clean(segment_b.get("Signal_diagnosis")),
    }
    diagnoses.discard("")
    if "Missing_or_unusable_alignment" in diagnoses or "Insufficient_recombinant_parent_overlap" in diagnoses:
        return "Insufficient_alignment_information"
    if "Low_phylogenetic_signal" in diagnoses:
        return "At_least_one_low_signal_segment"
    if "Closer_competing_sequence" in diagnoses:
        return "Possible_parent_mismatch_or_competing_lineage"
    if "Parent_indistinguishable_within_lineage" in diagnoses:
        return "Parental_lineage_resolved_but_exact_parent_ambiguous"
    if "Adequate_signal_but_unresolved_or_conflicting_branch" in diagnoses:
        return "Adequate_signal_but_topological_conflict"
    if diagnoses == {"Supported_segment"}:
        return "Both_segments_supported"
    return "Mixed_or_unclassified_signal"


def plot_diagnostics(
    segment_rows: Sequence[Mapping[str, object]],
    output_prefix: Path,
    formats: Sequence[str],
    minimum_ufboot: float,
) -> None:
    diagnoses = sorted({clean(row.get("Signal_diagnosis")) for row in segment_rows})
    palette = plt.get_cmap("tab10")
    colors = {diagnosis: palette(index % 10) for index, diagnosis in enumerate(diagnoses)}
    figure, axes = plt.subplots(1, 2, figsize=(14, 6))

    for diagnosis in diagnoses:
        selected = [row for row in segment_rows if clean(row.get("Signal_diagnosis")) == diagnosis]
        x = [to_float(row.get("Parsimony_informative_sites")) or 0.0 for row in selected]
        y = [to_float(row.get("UFBoot")) or 0.0 for row in selected]
        sizes = [
            25.0 + 8.0 * (to_int(row.get("Supported_bipartition_size")) or 1)
            for row in selected
        ]
        axes[0].scatter(x, y, s=sizes, alpha=0.75, label=diagnosis, color=colors[diagnosis])
    axes[0].axhline(minimum_ufboot, color="#d62728", linestyle="--", linewidth=1.2)
    axes[0].set_xlabel("Parsimony-informative sites")
    axes[0].set_ylabel("UFBoot (%)")
    axes[0].set_ylim(0, 103)
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=7, loc="best")

    counts = Counter(clean(row.get("Signal_diagnosis")) for row in segment_rows)
    labels = [item[0] for item in counts.most_common()]
    values = [counts[label] for label in labels]
    axes[1].barh(labels[::-1], values[::-1], color=[colors[label] for label in labels[::-1]])
    axes[1].set_xlabel("Number of segments")
    axes[1].set_title("Operational signal diagnosis")
    axes[1].grid(axis="x", alpha=0.2)

    figure.suptitle("Phylogenetic signal in insufficiently supported recombinant segments")
    figure.tight_layout()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        figure.savefig(f"{output_prefix}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_formats(value: str) -> List[str]:
    values = [item.strip().lower() for item in value.split(",") if item.strip()]
    invalid = sorted(set(values) - {"png", "pdf", "svg"})
    if invalid:
        raise argparse.ArgumentTypeError(f"Unsupported plot formats: {', '.join(invalid)}")
    return values


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quantify phylogenetic signal in insufficiently supported recombinant segments.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tree-summary", required=True, type=Path)
    parser.add_argument("--strict-calls", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--target-classes",
        default=DEFAULT_TARGET_CLASS,
        help="Comma-separated strict classes, or ALL",
    )
    parser.add_argument("--minimum-site-occupancy", type=float, default=0.70)
    parser.add_argument("--minimum-informative-sites", type=int, default=20)
    parser.add_argument("--minimum-informative-percent", type=float, default=1.0)
    parser.add_argument("--minimum-comparable-sites", type=int, default=50)
    parser.add_argument("--equivalent-identity-margin", type=float, default=0.5)
    parser.add_argument("--minimum-alrt", type=float, default=80.0)
    parser.add_argument("--minimum-ufboot", type=float, default=95.0)
    parser.add_argument("--formats", type=parse_formats, default=parse_formats("png,pdf"))
    args = parser.parse_args(argv)
    if not 0 < args.minimum_site_occupancy <= 1:
        parser.error("--minimum-site-occupancy must be in (0,1]")
    if args.minimum_informative_sites < 0 or args.minimum_comparable_sites < 0:
        parser.error("site count thresholds must be non-negative")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        _, tree_rows = read_tsv(args.tree_summary)
        _, call_rows = read_tsv(args.strict_calls)
        requested = {item.strip() for item in args.target_classes.split(",") if item.strip()}
        if "ALL" not in requested:
            call_rows = [row for row in call_rows if clean(row.get("Strict_validation_class")) in requested]
        if not call_rows:
            raise ValueError("No strict validation call remains after class filtering")

        tree_index = {key(row, include_segment=True): row for row in tree_rows}
        segment_rows: List[Dict[str, object]] = []
        segment_by_call: Dict[Tuple[str, ...], Dict[str, Dict[str, object]]] = defaultdict(dict)

        for call in call_rows:
            for segment in ("A_between_breakpoints", "B_complement"):
                tree_key = (*key(call), segment)
                tree_row = tree_index.get(tree_key)
                if tree_row is None:
                    raise ValueError(f"No tree summary row for {tree_key}")
                treefile = Path(clean(tree_row.get("Treefile")))
                fasta = treefile.with_suffix(".fasta")
                if not fasta.is_file():
                    raise FileNotFoundError(f"Segment FASTA not found: {fasta}")
                sequences = read_fasta(fasta)
                recombinant = first_token(call.get("Recombinant_accession"))
                parent = first_token(tree_row.get("Nearest_selected_parent"))
                metrics: Dict[str, object] = {
                    "Evidence_tier": call.get("Evidence_tier", ""),
                    "Recombinant_accession": recombinant,
                    "3SEQ_Event_ID": call.get("3SEQ_Event_ID", ""),
                    "RDP_Event_ID": call.get("RDP_Event_ID", ""),
                    "Breakpoint_scheme": call.get("Breakpoint_scheme", ""),
                    "Segment": segment,
                    "Strict_validation_class": call.get("Strict_validation_class", ""),
                    "Nearest_selected_parent": parent,
                    "Nearest_parent_roles": tree_row.get("Nearest_parent_roles", ""),
                    "SH_aLRT": tree_row.get("Nearest_parent_SH_aLRT", ""),
                    "UFBoot": tree_row.get("Nearest_parent_UFBoot", ""),
                    "Supported_bipartition_size": tree_row.get(
                        "Nearest_parent_supported_bipartition_size", ""
                    ),
                    "Tree_number_taxa": tree_row.get("Tree_number_taxa", ""),
                    "Segment_fasta": str(fasta),
                }
                alrt = to_float(metrics["SH_aLRT"])
                ufboot = to_float(metrics["UFBoot"])
                metrics["Branch_support_pass"] = support_passes(
                    alrt, ufboot, args.minimum_alrt, args.minimum_ufboot
                )
                metrics.update(alignment_site_metrics(sequences, args.minimum_site_occupancy))
                metrics.update(
                    parent_competition_metrics(
                        sequences,
                        recombinant,
                        parent,
                        args.equivalent_identity_margin,
                    )
                )
                metrics["Signal_diagnosis"] = diagnose_segment(
                    metrics,
                    args.minimum_informative_sites,
                    args.minimum_informative_percent,
                    args.minimum_comparable_sites,
                    args.equivalent_identity_margin,
                )
                segment_rows.append(metrics)
                segment_by_call[key(call)][segment] = metrics

        call_summaries: List[Dict[str, object]] = []
        for call in call_rows:
            segments = segment_by_call[key(call)]
            segment_a = segments.get("A_between_breakpoints", {})
            segment_b = segments.get("B_complement", {})
            call_summaries.append(
                {
                    "Evidence_tier": call.get("Evidence_tier", ""),
                    "Recombinant_accession": call.get("Recombinant_accession", ""),
                    "3SEQ_Event_ID": call.get("3SEQ_Event_ID", ""),
                    "RDP_Event_ID": call.get("RDP_Event_ID", ""),
                    "Breakpoint_scheme": call.get("Breakpoint_scheme", ""),
                    "Strict_validation_class": call.get("Strict_validation_class", ""),
                    "Role_switch_evidence": call.get("Role_switch_evidence", ""),
                    "Segment_A_signal_diagnosis": segment_a.get("Signal_diagnosis", ""),
                    "Segment_A_informative_sites": segment_a.get(
                        "Parsimony_informative_sites", ""
                    ),
                    "Segment_A_informative_percent": segment_a.get(
                        "Informative_percent_of_eligible", ""
                    ),
                    "Segment_A_parent_identity_percent": segment_a.get(
                        "Recombinant_parent_identity_percent", ""
                    ),
                    "Segment_A_parent_rank": segment_a.get("Parent_identity_rank", ""),
                    "Segment_B_signal_diagnosis": segment_b.get("Signal_diagnosis", ""),
                    "Segment_B_informative_sites": segment_b.get(
                        "Parsimony_informative_sites", ""
                    ),
                    "Segment_B_informative_percent": segment_b.get(
                        "Informative_percent_of_eligible", ""
                    ),
                    "Segment_B_parent_identity_percent": segment_b.get(
                        "Recombinant_parent_identity_percent", ""
                    ),
                    "Segment_B_parent_rank": segment_b.get("Parent_identity_rank", ""),
                    "Combined_signal_diagnosis": combined_diagnosis(segment_a, segment_b),
                }
            )

        args.output_dir.mkdir(parents=True, exist_ok=True)
        segment_output = args.output_dir / "segment_phylogenetic_signal.tsv"
        call_output = args.output_dir / "event_phylogenetic_signal_summary.tsv"
        write_tsv(segment_output, segment_rows, SEGMENT_FIELDS)
        write_tsv(call_output, call_summaries, CALL_FIELDS)
        plot_diagnostics(
            segment_rows,
            args.output_dir / "phylogenetic_signal_diagnostics",
            args.formats,
            args.minimum_ufboot,
        )

        segment_counts = Counter(clean(row["Signal_diagnosis"]) for row in segment_rows)
        call_counts = Counter(clean(row["Combined_signal_diagnosis"]) for row in call_summaries)
        print(f"Targeted event-schemes : {len(call_rows)}")
        print(f"Segments analyzed      : {len(segment_rows)}")
        print("\nSegment diagnoses:")
        for diagnosis, count in segment_counts.most_common():
            print(f"  {diagnosis:55s} {count:>3}")
        print("\nCombined diagnoses:")
        for diagnosis, count in call_counts.most_common():
            print(f"  {diagnosis:55s} {count:>3}")
        print(f"\nSegment table : {segment_output}")
        print(f"Event summary : {call_output}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
