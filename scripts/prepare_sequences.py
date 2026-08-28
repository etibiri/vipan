#!/usr/bin/env python3
"""Normalize, quality-filter and optionally rotate circular viral genomes.

The script deliberately writes an audit manifest for every input record.  It is
intended to be the first deterministic step of a Snakemake workflow, before any
multiple-sequence alignment is attempted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

from fasta_utils import FastaRecord, read_fasta, reverse_complement, write_fasta


VALID_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class MotifHit:
    orientation: str
    start: int
    mismatches: int


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def clean_sequence(value: str) -> str:
    return "".join(value.split()).upper().replace("U", "T")


def accession_from_record(record: FastaRecord) -> str:
    accession = record.id.strip()
    if not accession:
        raise ValueError("a FASTA record has an empty identifier")
    return accession


def circular_motif_hits(sequence: str, motif: str) -> Iterable[Tuple[int, int]]:
    """Yield ``(mismatches, start)`` for every circular motif placement."""

    if not sequence or len(sequence) < len(motif):
        return
    extended = sequence + sequence[: len(motif) - 1]
    for start in range(len(sequence)):
        window = extended[start : start + len(motif)]
        mismatches = sum(a != b for a, b in zip(window, motif))
        yield mismatches, start


def best_motif_hit(sequence: str, motif: str) -> MotifHit:
    candidates = []
    for orientation, oriented in (
        ("forward", sequence),
        ("reverse_complement", reverse_complement(sequence)),
    ):
        for mismatches, start in circular_motif_hits(oriented, motif):
            candidates.append((mismatches, 0 if orientation == "forward" else 1, start, orientation))
    if not candidates:
        return MotifHit("not_found", -1, len(motif) + 1)
    mismatches, _, start, orientation = min(candidates)
    return MotifHit(orientation, start, mismatches)


def orient_and_rotate(sequence: str, hit: MotifHit) -> str:
    if hit.orientation == "reverse_complement":
        sequence = reverse_complement(sequence)
    if hit.start < 0:
        return sequence
    return sequence[hit.start :] + sequence[: hit.start]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize and quality-filter circular viral FASTA sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--excluded-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--min-length", type=int, default=2400)
    parser.add_argument("--max-ambiguous-fraction", type=float, default=0.05)
    parser.add_argument("--orient-by-motif", action="store_true")
    parser.add_argument("--circular-motif", default="TAATATTAC")
    parser.add_argument("--max-motif-mismatches", type=int, default=1)
    parser.add_argument("--require-motif", action="store_true")
    parser.add_argument("--deduplicate-exact", action="store_true")
    args = parser.parse_args(argv)
    if args.min_length < 1:
        parser.error("--min-length must be positive")
    if not 0 <= args.max_ambiguous_fraction <= 1:
        parser.error("--max-ambiguous-fraction must be in [0,1]")
    if args.max_motif_mismatches < 0:
        parser.error("--max-motif-mismatches must be non-negative")
    args.circular_motif = clean_sequence(args.circular_motif)
    if not args.circular_motif or set(args.circular_motif) - VALID_BASES:
        parser.error("--circular-motif must contain only A, C, G and T")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        records = read_fasta(args.input)
        if not records:
            raise ValueError(f"no FASTA sequence found in {args.input}")

        seen_accessions = set()
        seen_sequences = {}
        kept = []
        excluded = []
        rows = []

        for source_order, record in enumerate(records, start=1):
            accession = accession_from_record(record)
            if accession in seen_accessions:
                raise ValueError(f"duplicate FASTA identifier: {accession}")
            seen_accessions.add(accession)

            sequence = clean_sequence(record.sequence)
            raw_length = len(sequence)
            ambiguous = sum(base not in VALID_BASES for base in sequence)
            ambiguous_fraction = ambiguous / raw_length if raw_length else 1.0
            hit = best_motif_hit(sequence, args.circular_motif)
            motif_pass = hit.mismatches <= args.max_motif_mismatches

            oriented = (
                orient_and_rotate(sequence, hit)
                if args.orient_by_motif and motif_pass
                else sequence
            )
            digest = sha256_text(oriented)

            status = "kept"
            reason = "passed_all_filters"
            representative = accession
            if raw_length < args.min_length:
                status = "excluded"
                reason = "sequence_too_short"
            elif ambiguous_fraction > args.max_ambiguous_fraction:
                status = "excluded"
                reason = "too_many_ambiguous_bases"
            elif args.require_motif and not motif_pass:
                status = "excluded"
                reason = "circular_origin_motif_not_found"
            elif args.deduplicate_exact and digest in seen_sequences:
                status = "excluded"
                reason = "exact_duplicate"
                representative = seen_sequences[digest]

            output_record = FastaRecord(accession, accession, oriented)
            if status == "kept":
                kept.append(output_record)
                seen_sequences.setdefault(digest, accession)
            else:
                excluded.append(output_record)

            rows.append(
                {
                    "source_order": source_order,
                    "accession": accession,
                    "original_description": record.description,
                    "raw_length": raw_length,
                    "ambiguous_bases": ambiguous,
                    "ambiguous_fraction": f"{ambiguous_fraction:.8f}",
                    "motif": args.circular_motif,
                    "motif_orientation": hit.orientation,
                    "motif_start_zero_based": hit.start if hit.start >= 0 else "",
                    "motif_mismatches": hit.mismatches,
                    "sequence_sha256": digest,
                    "status": status,
                    "reason": reason,
                    "representative_accession": representative,
                }
            )

        if len(kept) < 3:
            raise ValueError(f"only {len(kept)} sequences passed preprocessing; at least 3 are required")

        for path in (args.output_fasta, args.excluded_fasta, args.manifest):
            path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(kept, args.output_fasta)
        write_fasta(excluded, args.excluded_fasta)
        with args.manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        print(f"Input sequences    : {len(records)}")
        print(f"Kept sequences     : {len(kept)}")
        print(f"Excluded sequences : {len(excluded)}")
        print(f"Manifest           : {args.manifest}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
