#!/usr/bin/env python3
"""Apply strict, reproducible classes to segment-tree validation calls."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def as_float(value: object) -> Optional[float]:
    try:
        return float(clean(value))
    except ValueError:
        return None


def as_int(value: object) -> Optional[int]:
    number = as_float(value)
    return int(number) if number is not None else None


def as_bool(value: object) -> bool:
    return clean(value).lower() in {"true", "t", "yes", "y", "1"}


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def event_key(row: Mapping[str, object]) -> Tuple[str, str, str, str]:
    return (
        clean(row.get("Recombinant_accession")),
        clean(row.get("3SEQ_Event_ID")),
        clean(row.get("RDP_Event_ID")),
        clean(row.get("Breakpoint_scheme")),
    )


def tree_key(row: Mapping[str, object]) -> Tuple[str, str, str, str, str]:
    return (*event_key(row), clean(row.get("Segment")))


def support_passes(row: Mapping[str, object], minimum_alrt: float, minimum_ufboot: float) -> bool:
    alrt = as_float(row.get("Nearest_parent_SH_aLRT"))
    ufboot = as_float(row.get("Nearest_parent_UFBoot"))
    return alrt is not None and ufboot is not None and alrt >= minimum_alrt and ufboot >= minimum_ufboot


CLASS_RANK = {
    "A_local_supported_parent_switch": 0,
    "B_intermediate_supported_parent_switch": 1,
    "C_broad_lineage_supported_switch": 2,
    "D_expected_switch_insufficient_branch_support": 3,
    "E_parent_switch_without_expected_roles": 4,
    "F_no_parent_switch_or_unavailable": 5,
}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Classify recombination validation calls using branch support and bipartition size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--calls", required=True, type=Path)
    parser.add_argument("--tree-summary", required=True, type=Path)
    parser.add_argument("--output-calls", required=True, type=Path)
    parser.add_argument("--output-events", required=True, type=Path)
    parser.add_argument("--minimum-alrt", type=float, default=80.0)
    parser.add_argument("--minimum-ufboot", type=float, default=95.0)
    parser.add_argument("--local-bipartition-maximum", type=int, default=4)
    parser.add_argument("--intermediate-bipartition-maximum", type=int, default=6)
    args = parser.parse_args(argv)
    if args.local_bipartition_maximum < 2:
        parser.error("--local-bipartition-maximum must be at least 2")
    if args.intermediate_bipartition_maximum < args.local_bipartition_maximum:
        parser.error("intermediate maximum must be >= local maximum")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        call_fields, calls = read_tsv(args.calls)
        _, trees = read_tsv(args.tree_summary)
        tree_index = {tree_key(row): row for row in trees}
        detailed = []

        for call in calls:
            key = event_key(call)
            segment_a = tree_index.get((*key, "A_between_breakpoints"), {})
            segment_b = tree_index.get((*key, "B_complement"), {})
            a_pass = support_passes(segment_a, args.minimum_alrt, args.minimum_ufboot)
            b_pass = support_passes(segment_b, args.minimum_alrt, args.minimum_ufboot)
            size_a = as_int(segment_a.get("Nearest_parent_supported_bipartition_size"))
            size_b = as_int(segment_b.get("Nearest_parent_supported_bipartition_size"))
            expected = as_bool(call.get("Expected_role_switch"))
            different = as_bool(call.get("Different_nearest_parents"))

            if expected and a_pass and b_pass and size_a is not None and size_b is not None:
                maximum_size = max(size_a, size_b)
                if maximum_size <= args.local_bipartition_maximum:
                    strict_class = "A_local_supported_parent_switch"
                elif maximum_size <= args.intermediate_bipartition_maximum:
                    strict_class = "B_intermediate_supported_parent_switch"
                else:
                    strict_class = "C_broad_lineage_supported_switch"
            elif expected:
                strict_class = "D_expected_switch_insufficient_branch_support"
            elif different:
                strict_class = "E_parent_switch_without_expected_roles"
            else:
                strict_class = "F_no_parent_switch_or_unavailable"

            detailed.append(
                {
                    **call,
                    "Segment_A_supported_bipartition_size": "" if size_a is None else size_a,
                    "Segment_B_supported_bipartition_size": "" if size_b is None else size_b,
                    "Segment_A_branch_support_pass": a_pass,
                    "Segment_B_branch_support_pass": b_pass,
                    "Strict_validation_class": strict_class,
                    "Strict_minimum_SH_aLRT": args.minimum_alrt,
                    "Strict_minimum_UFBoot": args.minimum_ufboot,
                    "Strict_local_bipartition_maximum": args.local_bipartition_maximum,
                    "Strict_intermediate_bipartition_maximum": args.intermediate_bipartition_maximum,
                }
            )

        grouped: Dict[str, List[Mapping[str, object]]] = defaultdict(list)
        for row in detailed:
            grouped[clean(row.get("Recombinant_accession"))].append(row)
        event_rows = []
        for recombinant, rows in sorted(grouped.items()):
            best = min(rows, key=lambda row: (CLASS_RANK[clean(row["Strict_validation_class"])], clean(row.get("Breakpoint_scheme"))))
            event_rows.append(
                {
                    "Recombinant_accession": recombinant,
                    "Evidence_tier": best.get("Evidence_tier", ""),
                    "3SEQ_Event_ID": best.get("3SEQ_Event_ID", ""),
                    "RDP_Event_ID": best.get("RDP_Event_ID", ""),
                    "Best_breakpoint_scheme": best.get("Breakpoint_scheme", ""),
                    "Best_strict_validation_class": best.get("Strict_validation_class", ""),
                    "Number_schemes_evaluated": len(rows),
                    "All_scheme_classes": ";".join(
                        f"{clean(row.get('Breakpoint_scheme'))}:{clean(row.get('Strict_validation_class'))}"
                        for row in sorted(rows, key=lambda item: clean(item.get("Breakpoint_scheme")))
                    ),
                }
            )

        extra_fields = [
            "Segment_A_supported_bipartition_size",
            "Segment_B_supported_bipartition_size",
            "Segment_A_branch_support_pass",
            "Segment_B_branch_support_pass",
            "Strict_validation_class",
            "Strict_minimum_SH_aLRT",
            "Strict_minimum_UFBoot",
            "Strict_local_bipartition_maximum",
            "Strict_intermediate_bipartition_maximum",
        ]
        write_tsv(args.output_calls, detailed, [*call_fields, *extra_fields])
        write_tsv(args.output_events, event_rows, list(event_rows[0]))

        scheme_counts = Counter(clean(row["Strict_validation_class"]) for row in detailed)
        event_counts = Counter(clean(row["Best_strict_validation_class"]) for row in event_rows)
        print("Strict classes by breakpoint scheme:")
        for label, count in scheme_counts.most_common():
            print(f"  {label:55s} {count:>4}")
        print("Strict classes by recombinant:")
        for label, count in event_counts.most_common():
            print(f"  {label:55s} {count:>4}")
        print(f"Detailed calls : {args.output_calls}")
        print(f"Event summary  : {args.output_events}")
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
