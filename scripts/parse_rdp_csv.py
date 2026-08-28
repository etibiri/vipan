#!/usr/bin/env python3
"""Parse the hierarchical RDP4 CSV report into audited event tables.

The parser supports continuation rows used by RDP for events containing several
recombinants or parental candidates.  Robustness decisions are explicit command
line parameters and are copied to every output row.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from fasta_utils import read_fasta


ACCESSION_RE = re.compile(
    r"(?<![A-Za-z0-9_])(?:[A-Z]{1,5}_\d{4,10}|[A-Z]{1,5}\d{5,10})"
    r"(?:\.\d+)?(?![A-Za-z0-9_])",
    flags=re.IGNORECASE,
)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def parse_int(value: object) -> Optional[int]:
    match = re.search(r"\d+", clean(value))
    return int(match.group(0)) if match else None


def parse_pvalue(value: object) -> Optional[float]:
    text = clean(value).replace(",", "")
    if not text or text.upper() in {"NS", "NA", "N/A", "-"}:
        return None
    text = text.lstrip("~*^")
    try:
        return float(text)
    except ValueError:
        return None


def extract_accessions(value: object) -> List[str]:
    text = clean(value).strip("\"'")
    if not text or text.lower() in {"unknown", "none", "na", "n/a", "-"}:
        return []
    found = [match.group(0).upper() for match in ACCESSION_RE.finditer(text)]
    if found:
        return list(dict.fromkeys(found))
    tokens = []
    for part in re.split(r"[;,]+", text):
        token = clean(part).lstrip("~*^").split()[0] if clean(part) else ""
        if token and token.lower() not in {"unknown", "none", "na", "n/a", "-"}:
            tokens.append(token.upper())
    return list(dict.fromkeys(tokens))


def unique(values: Iterable[str]) -> List[str]:
    return list(dict.fromkeys(value for value in values if value))


def canonical_method(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", clean(value).lower())


def alignment_info(path: Path) -> Tuple[int, set]:
    records = read_fasta(path)
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise ValueError("--alignment must be an aligned FASTA")
    return lengths.pop(), {record.id.split()[0].upper() for record in records}


def read_report(path: Path) -> Tuple[List[str], List[List[str]], List[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))
    header_index = None
    for index, row in enumerate(rows):
        if row and clean(row[0]).lower().startswith("recombination event number"):
            header_index = index
            break
    if header_index is None:
        raise ValueError("RDP event header was not found")
    header = rows[header_index]
    methods = [clean(value) for value in header[11:] if clean(value)]
    data_rows = []
    notes = []
    for row in rows[header_index + 1 :]:
        if row and clean(row[0]).lower() == "table key:":
            break
        if not row or not any(clean(value) for value in row):
            continue
        if len(row) == 1:
            notes.append(clean(row[0]))
            continue
        data_rows.append(row)
    return methods, data_rows, notes


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize and robustly filter an RDP4 CSV event report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--output-all", required=True, type=Path)
    parser.add_argument("--output-robust", required=True, type=Path)
    parser.add_argument("--min-methods", type=int, default=3)
    parser.add_argument("--min-tested-methods", type=int, default=7)
    parser.add_argument("--min-support-fraction", type=float, default=0.50)
    parser.add_argument("--max-missing-methods", type=int, default=2)
    parser.add_argument("--method-pvalue-cutoff", type=float, default=0.05)
    parser.add_argument("--min-fragment-length", type=int, default=100)
    parser.add_argument("--include-internal-3seq", action="store_true")
    parser.add_argument("--require-known-parents", action="store_true")
    parser.add_argument("--allow-caret", action="store_true")
    parser.add_argument("--allow-tilde", action="store_true")
    args = parser.parse_args(argv)
    if args.min_methods < 1 or args.min_tested_methods < 1:
        parser.error("method-count thresholds must be positive")
    if not 0 <= args.min_support_fraction <= 1:
        parser.error("--min-support-fraction must be in [0,1]")
    if not 0 < args.method_pvalue_cutoff <= 1:
        parser.error("--method-pvalue-cutoff must be in (0,1]")
    return args


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fragment_length(begin: Optional[int], end: Optional[int], length: int) -> Optional[int]:
    if begin is None or end is None:
        return None
    if begin <= end:
        return end - begin + 1
    return length - begin + 1 + end


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        alignment_length, alignment_accessions = alignment_info(args.alignment)
        method_names, raw_rows, _ = read_report(args.input)
        if not method_names:
            raise ValueError("no RDP detection-method columns were found")

        grouped: Dict[str, Dict[str, object]] = {}
        current_event = ""
        for source_line, raw in enumerate(raw_rows, start=1):
            row = list(raw) + [""] * max(0, 11 + len(method_names) - len(raw))
            event = clean(row[0]) or current_event
            if not event:
                continue
            current_event = event
            has_coordinates = any(clean(row[index]) for index in range(2, 8))
            has_methods = any(clean(value) for value in row[11 : 11 + len(method_names)])
            if event not in grouped:
                grouped[event] = {
                    "source_line": source_line,
                    "rdp_file_number": clean(row[1]),
                    "alignment_begin": clean(row[2]),
                    "alignment_end": clean(row[3]),
                    "recombinant_begin": clean(row[4]),
                    "recombinant_end": clean(row[5]),
                    "reference_begin": clean(row[6]),
                    "reference_end": clean(row[7]),
                    "recombinant_raw": [],
                    "minor_raw": [],
                    "major_raw": [],
                    "method_raw": {name: "" for name in method_names},
                }
            target = grouped[event]
            if has_coordinates or has_methods:
                for key, index in (
                    ("rdp_file_number", 1),
                    ("alignment_begin", 2),
                    ("alignment_end", 3),
                    ("recombinant_begin", 4),
                    ("recombinant_end", 5),
                    ("reference_begin", 6),
                    ("reference_end", 7),
                ):
                    if clean(row[index]):
                        target[key] = clean(row[index])
                for name, value in zip(method_names, row[11 : 11 + len(method_names)]):
                    if clean(value):
                        target["method_raw"][name] = clean(value)
            for key, index in (("recombinant_raw", 8), ("minor_raw", 9), ("major_raw", 10)):
                if clean(row[index]):
                    target[key].append(clean(row[index]))

        output_rows = []
        excluded_method_keys = set() if args.include_internal_3seq else {"3seq", "threeseq"}
        for event, item in sorted(grouped.items(), key=lambda pair: (parse_int(pair[0]) or math.inf, pair[0])):
            recombinants = unique(
                accession
                for value in item["recombinant_raw"]
                for accession in extract_accessions(value)
            )
            minor = unique(
                accession for value in item["minor_raw"] for accession in extract_accessions(value)
            )
            major = unique(
                accession for value in item["major_raw"] for accession in extract_accessions(value)
            )
            method_columns: Dict[str, object] = {}
            significant = nonsignificant = tested = not_reported = 0
            pvalues = []
            for name in method_names:
                key = canonical_method(name)
                raw_value = clean(item["method_raw"].get(name, ""))
                pvalue = parse_pvalue(raw_value)
                if not raw_value:
                    status = "not_reported"
                elif raw_value.upper() == "NS":
                    status = "nonsignificant"
                elif pvalue is None:
                    status = "unparsed"
                elif pvalue <= args.method_pvalue_cutoff:
                    status = "significant"
                else:
                    status = "nonsignificant"
                method_columns[f"{key}_raw"] = raw_value
                method_columns[f"{key}_status"] = status
                method_columns[f"{key}_pvalue"] = "" if pvalue is None else f"{pvalue:.12g}"
                if key in excluded_method_keys:
                    continue
                if status == "not_reported":
                    not_reported += 1
                elif status in {"significant", "nonsignificant"}:
                    tested += 1
                    if status == "significant":
                        significant += 1
                        if pvalue is not None:
                            pvalues.append(pvalue)
                    else:
                        nonsignificant += 1

            begin = parse_int(item["alignment_begin"])
            end = parse_int(item["alignment_end"])
            length = fragment_length(begin, end, alignment_length)
            raw_text = " ".join(
                [*item["recombinant_raw"], *item["minor_raw"], *item["major_raw"]]
            )
            contains_caret = "^" in raw_text
            contains_tilde = "~" in raw_text
            unknown_minor = any("unknown" in value.lower() for value in item["minor_raw"])
            unknown_major = any("unknown" in value.lower() for value in item["major_raw"])
            all_recombinants_in_alignment = bool(recombinants) and all(
                accession.upper() in alignment_accessions for accession in recombinants
            )
            known_parents = [*minor, *major]
            all_known_parents_in_alignment = bool(known_parents) and all(
                accession.upper() in alignment_accessions for accession in known_parents
            )
            support_fraction = significant / tested if tested else 0.0
            failed = []
            if significant < args.min_methods:
                failed.append("min_methods")
            if tested < args.min_tested_methods:
                failed.append("min_tested_methods")
            if support_fraction < args.min_support_fraction:
                failed.append("support_fraction")
            if not_reported > args.max_missing_methods:
                failed.append("max_missing_methods")
            if length is None or length < args.min_fragment_length:
                failed.append("min_fragment_length")
            if not all_recombinants_in_alignment:
                failed.append("recombinant_absent_from_alignment")
            if args.require_known_parents and (unknown_minor or unknown_major or not minor or not major):
                failed.append("known_parents_required")
            if not args.allow_caret and contains_caret:
                failed.append("caret_marker")
            if not args.allow_tilde and contains_tilde:
                failed.append("alternative_process_marker")

            output_rows.append(
                {
                    "event": event,
                    "rdp_file_number": item["rdp_file_number"],
                    "alignment_length": alignment_length,
                    "alignment_begin": "" if begin is None else begin,
                    "alignment_end": "" if end is None else end,
                    "crosses_alignment_origin": bool(begin is not None and end is not None and begin > end),
                    "recombinant_begin": item["recombinant_begin"],
                    "recombinant_end": item["recombinant_end"],
                    "reference_begin": item["reference_begin"],
                    "reference_end": item["reference_end"],
                    "reference_fragment_length": "" if length is None else length,
                    "number_recombinants": len(recombinants),
                    "recombinant_accessions": ";".join(recombinants),
                    "recombinant_sequences_raw": ";".join(item["recombinant_raw"]),
                    "number_minor_parent_candidates": len(minor),
                    "minor_parent_accessions": ";".join(minor),
                    "minor_parent_candidates_raw": ";".join(item["minor_raw"]),
                    "number_major_parent_candidates": len(major),
                    "major_parent_accessions": ";".join(major),
                    "major_parent_candidates_raw": ";".join(item["major_raw"]),
                    "all_recombinants_in_alignment": all_recombinants_in_alignment,
                    "all_known_parents_in_alignment": all_known_parents_in_alignment,
                    **method_columns,
                    "number_significant_methods": significant,
                    "number_nonsignificant_methods": nonsignificant,
                    "number_tested_methods": tested,
                    "number_not_reported_methods": not_reported,
                    "support_fraction": f"{support_fraction:.8f}",
                    "minimum_pvalue": "" if not pvalues else f"{min(pvalues):.12g}",
                    "contains_recombinant_caret_marker": contains_caret,
                    "contains_alternative_process_marker": contains_tilde,
                    "contains_unknown_minor_parent": unknown_minor,
                    "contains_unknown_major_parent": unknown_major,
                    "criterion_min_methods": args.min_methods,
                    "criterion_min_tested_methods": args.min_tested_methods,
                    "criterion_support_fraction": args.min_support_fraction,
                    "criterion_max_missing_methods": args.max_missing_methods,
                    "criterion_method_pvalue_cutoff": args.method_pvalue_cutoff,
                    "criterion_min_fragment_length": args.min_fragment_length,
                    "criterion_internal_3seq_excluded": not args.include_internal_3seq,
                    "criterion_require_known_parents": args.require_known_parents,
                    "criterion_allow_caret": args.allow_caret,
                    "criterion_allow_tilde": args.allow_tilde,
                    "robust_event": not failed,
                    "failed_robustness_criteria": ";".join(failed),
                }
            )

        if not output_rows:
            raise ValueError("no RDP events were parsed")
        robust_rows = [row for row in output_rows if row["robust_event"]]
        fields = list(output_rows[0])
        write_tsv(args.output_all, output_rows, fields)
        write_tsv(args.output_robust, robust_rows, fields)
        print(f"RDP events parsed : {len(output_rows)}")
        print(f"Robust events     : {len(robust_rows)}")
        print(f"All events        : {args.output_all}")
        print(f"Robust events     : {args.output_robust}")
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
