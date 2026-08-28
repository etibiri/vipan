#!/usr/bin/env python3
"""Create publication-ready QC summaries for a nucleotide alignment."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from fasta_utils import read_fasta


VALID_BASES = frozenset("ACGT")
GAPS = frozenset("-.")


def pairwise_p_distance(first: str, second: str) -> Optional[float]:
    comparable = 0
    different = 0
    for base_a, base_b in zip(first, second):
        if base_a in VALID_BASES and base_b in VALID_BASES:
            comparable += 1
            different += base_a != base_b
    return different / comparable if comparable else None


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize an aligned nucleotide FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--per-sequence-tsv", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--minimum-site-occupancy", type=float, default=0.70)
    args = parser.parse_args(argv)
    if not 0 < args.minimum_site_occupancy <= 1:
        parser.error("--minimum-site-occupancy must be in (0,1]")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        records = read_fasta(args.alignment)
        if len(records) < 2:
            raise ValueError("at least two aligned sequences are required")
        lengths = {len(record.sequence) for record in records}
        if len(lengths) != 1:
            raise ValueError(f"FASTA is not aligned; lengths={sorted(lengths)[:10]}")
        alignment_length = lengths.pop()
        number_sequences = len(records)
        minimum_valid = max(2, math.ceil(args.minimum_site_occupancy * number_sequences))
        sequences = {record.id.split()[0]: record.sequence.upper() for record in records}
        if len(sequences) != number_sequences:
            raise ValueError("duplicate sequence identifiers are present")

        per_sequence = []
        for accession, sequence in sorted(sequences.items()):
            valid = sum(base in VALID_BASES for base in sequence)
            gaps = sum(base in GAPS for base in sequence)
            ambiguity = alignment_length - valid - gaps
            per_sequence.append(
                {
                    "accession": accession,
                    "alignment_length": alignment_length,
                    "valid_bases": valid,
                    "valid_fraction": f"{valid/alignment_length:.8f}",
                    "gap_bases": gaps,
                    "gap_fraction": f"{gaps/alignment_length:.8f}",
                    "ambiguous_bases": ambiguity,
                    "ambiguity_fraction": f"{ambiguity/alignment_length:.8f}",
                }
            )

        eligible = constant = variable = singleton = informative = 0
        occupancy_values: List[float] = []
        for position in range(alignment_length):
            bases = [sequence[position] for sequence in sequences.values() if sequence[position] in VALID_BASES]
            occupancy_values.append(len(bases) / number_sequences)
            if len(bases) < minimum_valid:
                continue
            eligible += 1
            counts = Counter(bases)
            if len(counts) == 1:
                constant += 1
            else:
                variable += 1
                if sum(count >= 2 for count in counts.values()) >= 2:
                    informative += 1
                else:
                    singleton += 1

        pairwise = []
        for first, second in itertools.combinations(sorted(sequences), 2):
            current = pairwise_p_distance(sequences[first], sequences[second])
            if current is not None:
                pairwise.append(current)

        metrics: Dict[str, object] = {
            "number_sequences": number_sequences,
            "alignment_length": alignment_length,
            "minimum_valid_sequences_per_site": minimum_valid,
            "eligible_sites": eligible,
            "excluded_low_occupancy_sites": alignment_length - eligible,
            "mean_site_occupancy_percent": 100 * statistics.mean(occupancy_values),
            "constant_sites": constant,
            "variable_sites": variable,
            "singleton_variable_sites": singleton,
            "parsimony_informative_sites": informative,
            "variable_percent_of_eligible": 100 * variable / eligible if eligible else 0,
            "informative_percent_of_eligible": 100 * informative / eligible if eligible else 0,
            "unique_aligned_haplotypes": len(set(sequences.values())),
            "pairwise_p_distance_median": statistics.median(pairwise) if pairwise else "",
            "pairwise_p_distance_minimum": min(pairwise) if pairwise else "",
            "pairwise_p_distance_maximum": max(pairwise) if pairwise else "",
        }

        for path in (args.summary_tsv, args.per_sequence_tsv, args.summary_json):
            path.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_tsv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Metric", "Value"])
            for key, value in metrics.items():
                writer.writerow([key, f"{value:.8f}" if isinstance(value, float) else value])
        with args.per_sequence_tsv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(per_sequence[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(per_sequence)
        args.summary_json.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Sequences                    : {number_sequences}")
        print(f"Alignment length             : {alignment_length}")
        print(f"Variable sites               : {variable}")
        print(f"Parsimony-informative sites  : {informative}")
        print(f"Summary                      : {args.summary_tsv}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
