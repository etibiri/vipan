#!/usr/bin/env python3
"""Compare recombination events reported by standalone 3SEQ and RDP5.

The comparison is performed at the event level, not only at the recombinant
accession level.  Events are first grouped by normalized recombinant accession.
Within each accession, a minimum-cost one-to-one assignment is used so that an
RDP event cannot confirm several 3SEQ events (or conversely).

The program supports:

* RDP5 robust-event TSV files produced by ``parse_rdp_report.py``;
* non-redundant 3SEQ event TSV files;
* native ``.3s.rec`` files, including rows that contain more breakpoint fields
  than the nominal 13-column header;
* multiple alternative breakpoint-range pairs reported by 3SEQ;
* linear or circular coordinates;
* direct or swapped breakpoint ordering;
* separate assessment of breakpoint and parental concordance.

Outputs are written as:

* ``<prefix>.event_pairs.tsv``: detailed one-to-one event comparisons;
* ``<prefix>.recombinant_summary.tsv``: per-recombinant summary;
* ``<prefix>.summary.tsv``: global metrics and concordance counts.

Only Python's standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


ACCESSION_RE = re.compile(
    r"(?<![A-Za-z0-9_])(?:[A-Z]{1,5}_\d{4,10}|[A-Z]{1,5}\d{5,10})"
    r"(?:\.\d+)?(?![A-Za-z0-9_])",
    flags=re.IGNORECASE,
)

BREAKPOINT_PAIR_RE = re.compile(
    r"(?P<a>\d+)\s*(?:-\s*(?P<b>\d+))?\s*&\s*"
    r"(?P<c>\d+)\s*(?:-\s*(?P<d>\d+))?"
)

SINGLE_INTERVAL_RE = re.compile(
    r"^\s*[~*^]?\s*(?P<a>\d+)\s*(?:-\s*(?P<b>\d+))?\s*$"
)


def eprint(*args: object) -> None:
    """Print a diagnostic message to stderr."""

    print(*args, file=sys.stderr)


def canonical_name(value: object) -> str:
    """Return a stable lowercase representation of a column name."""

    text = str(value or "").strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def clean_text(value: object) -> str:
    """Normalize spaces and non-breaking spaces without altering identifiers."""

    return re.sub(r"\s+", " ", str(value or "").replace("\u00a0", " ")).strip()


def get_first(row: Mapping[str, str], *names: str) -> str:
    """Return the first non-empty value for a list of canonical aliases."""

    for name in names:
        value = row.get(canonical_name(name), "")
        if clean_text(value):
            return clean_text(value)
    return ""


def parse_bool(value: object) -> Optional[bool]:
    """Parse common textual boolean values."""

    text = clean_text(value).lower()
    if text in {"true", "t", "yes", "y", "1", "oui"}:
        return True
    if text in {"false", "f", "no", "n", "0", "non"}:
        return False
    return None


def parse_float(value: object) -> Optional[float]:
    """Parse a float while accepting infinities and scientific notation."""

    text = clean_text(value)
    if not text or text.upper() in {"NA", "N/A", "NONE", "."}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_int(value: object) -> Optional[int]:
    """Extract the first integer from a possibly marked coordinate field."""

    match = re.search(r"\d+", clean_text(value))
    return int(match.group(0)) if match else None


def format_number(value: Optional[float], digits: int = 6) -> str:
    """Format an optional numeric value for TSV output."""

    if value is None:
        return ""
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return f"{value:.{digits}f}"


def extract_accessions(value: object) -> List[str]:
    """Extract normalized accession-like identifiers from a text field.

    RDP continuation rows may result in several accessions in one field.  If no
    accession pattern is found, delimited tokens are retained as a fallback so
    that non-INSDC sequence identifiers remain usable.
    """

    text = clean_text(value).strip("\"'")
    if not text or text.lower() in {"unknown", "none", "na", "n/a", "-"}:
        return []

    found = [m.group(0).upper() for m in ACCESSION_RE.finditer(text)]
    if found:
        return list(dict.fromkeys(found))

    tokens: List[str] = []
    for token in re.split(r"[;,]+", text):
        token = clean_text(token).lstrip("^~*").strip("\"'")
        if not token or token.lower() in {"unknown", "none", "na", "n/a", "-"}:
            continue
        token = token.split()[0]
        if "|" in token:
            pipe_tokens = [part for part in token.split("|") if part]
            token = pipe_tokens[-1] if pipe_tokens else token
        tokens.append(token.upper())
    return list(dict.fromkeys(tokens))


def first_accession(value: object) -> str:
    """Return the first normalized sequence identifier in a field."""

    accessions = extract_accessions(value)
    return accessions[0] if accessions else ""


def detect_delimiter(path: Path) -> str:
    """Detect tab- versus comma-delimited text from the first non-empty line."""

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.strip():
                return "\t" if line.count("\t") >= line.count(",") else ","
    raise ValueError(f"Empty input file: {path}")


def unique_headers(header: Sequence[str]) -> List[str]:
    """Canonicalize headers and disambiguate duplicate names."""

    counts: Counter[str] = Counter()
    result: List[str] = []
    for index, raw_name in enumerate(header, start=1):
        base = canonical_name(raw_name) or f"column_{index}"
        counts[base] += 1
        result.append(base if counts[base] == 1 else f"{base}_{counts[base]}")
    return result


def read_table(path: Path) -> Tuple[List[Dict[str, str]], List[str], str]:
    """Read a delimited table while preserving surplus columns.

    Native 3SEQ output can declare 13 header fields but contain 15 fields in a
    data row.  Surplus values are retained under ``__extra_1``,
    ``__extra_2``, etc., rather than being discarded.
    """

    delimiter = detect_delimiter(path)
    rows: List[Dict[str, str]] = []

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header_raw: Optional[List[str]] = None
        header: List[str] = []

        for values in reader:
            if not values or not any(clean_text(v) for v in values):
                continue
            if header_raw is None:
                header_raw = values
                header = unique_headers(values)
                continue

            row: Dict[str, str] = {}
            for index, name in enumerate(header):
                row[name] = values[index] if index < len(values) else ""
            for extra_index, value in enumerate(values[len(header) :], start=1):
                row[f"__extra_{extra_index}"] = value
            rows.append(row)

    if header_raw is None:
        raise ValueError(f"No header found in input file: {path}")
    return rows, header, delimiter


@dataclass(frozen=True)
class Interval:
    """A one-based inclusive breakpoint interval."""

    start: int
    end: int

    def label(self) -> str:
        return str(self.start) if self.start == self.end else f"{self.start}-{self.end}"


@dataclass(frozen=True)
class BreakpointAlternative:
    """One possible pair of breakpoint ranges reported by 3SEQ."""

    breakpoint_1: Interval
    breakpoint_2: Interval
    raw: str
    source: str

    def label(self) -> str:
        return f"{self.breakpoint_1.label()} & {self.breakpoint_2.label()}"


@dataclass
class ThreeSeqEvent:
    event_id: str
    recombinant: str
    parent_1: str
    parent_2: str
    alternatives: List[BreakpointAlternative]
    p_value: Optional[float]
    min_recombinant_length: Optional[int]
    source_row: int


@dataclass
class RDPEvent:
    event_id: str
    recombinant: str
    minor_parents: Tuple[str, ...]
    major_parents: Tuple[str, ...]
    breakpoint_1: Optional[int]
    breakpoint_2: Optional[int]
    crosses_origin: Optional[bool]
    p_value: Optional[float]
    source_row: int


@dataclass
class AlternativeComparison:
    alternative_index: int
    alternative: BreakpointAlternative
    orientation: str
    matched_rdp_for_bp1: Optional[int]
    matched_rdp_for_bp2: Optional[int]
    distance_1: Optional[float]
    distance_2: Optional[float]
    maximum_distance: Optional[float]
    sum_distance: Optional[float]
    partition_jaccard: Optional[float]


@dataclass
class EventComparison:
    recombinant: str
    threeseq: Optional[ThreeSeqEvent]
    rdp: Optional[RDPEvent]
    match_status: str
    concordance_class: str
    parent_concordance: str
    best: Optional[AlternativeComparison]
    assignment_cost: Optional[float]


def parse_single_interval(value: object) -> Optional[Interval]:
    """Parse an integer or ``start-end`` breakpoint interval."""

    match = SINGLE_INTERVAL_RE.match(clean_text(value))
    if not match:
        return None
    start = int(match.group("a"))
    end = int(match.group("b") or match.group("a"))
    return Interval(start, end)


def parse_breakpoint_pairs(text: object, source: str) -> List[BreakpointAlternative]:
    """Extract all ``range & range`` alternatives from a text value."""

    cleaned = clean_text(text)
    alternatives: List[BreakpointAlternative] = []
    for match in BREAKPOINT_PAIR_RE.finditer(cleaned):
        a = int(match.group("a"))
        b = int(match.group("b") or match.group("a"))
        c = int(match.group("c"))
        d = int(match.group("d") or match.group("c"))
        raw = clean_text(match.group(0))
        alternatives.append(
            BreakpointAlternative(Interval(a, b), Interval(c, d), raw, source)
        )
    return alternatives


def deduplicate_alternatives(
    alternatives: Iterable[BreakpointAlternative],
) -> List[BreakpointAlternative]:
    """Remove duplicate breakpoint pairs while retaining input order."""

    seen = set()
    result: List[BreakpointAlternative] = []
    for alternative in alternatives:
        key = (
            alternative.breakpoint_1.start,
            alternative.breakpoint_1.end,
            alternative.breakpoint_2.start,
            alternative.breakpoint_2.end,
        )
        if key not in seen:
            seen.add(key)
            result.append(alternative)
    return result


def extract_3seq_alternatives(row: Mapping[str, str]) -> List[BreakpointAlternative]:
    """Extract all 3SEQ breakpoint alternatives from normalized or native rows."""

    alternatives: List[BreakpointAlternative] = []

    # Prefer explicit raw/native breakpoint fields because they retain ranges.
    for name, value in row.items():
        if not clean_text(value):
            continue
        if (
            name.startswith("__extra_")
            or name in {"breakpoints", "breakpoint_range_raw", "breakpoint_ranges"}
            or ("breakpoint" in name and "raw" in name)
        ):
            alternatives.extend(parse_breakpoint_pairs(value, source=name))

    if alternatives:
        return deduplicate_alternatives(alternatives)

    # Fall back to representative breakpoint columns in a normalized event file.
    breakpoint_1 = parse_single_interval(
        get_first(row, "breakpoint_1", "breakpoint1", "bp1")
    )
    breakpoint_2 = parse_single_interval(
        get_first(row, "breakpoint_2", "breakpoint2", "bp2")
    )
    if breakpoint_1 and breakpoint_2:
        alternatives.append(
            BreakpointAlternative(
                breakpoint_1,
                breakpoint_2,
                f"{breakpoint_1.label()} & {breakpoint_2.label()}",
                "representative_breakpoints",
            )
        )

    # Cluster centers are only used when no representative/range pair is available.
    if not alternatives:
        cluster_1 = parse_single_interval(
            get_first(row, "breakpoint_cluster_1", "breakpoint_cluster1")
        )
        cluster_2 = parse_single_interval(
            get_first(row, "breakpoint_cluster_2", "breakpoint_cluster2")
        )
        if cluster_1 and cluster_2:
            alternatives.append(
                BreakpointAlternative(
                    cluster_1,
                    cluster_2,
                    f"{cluster_1.label()} & {cluster_2.label()}",
                    "breakpoint_clusters",
                )
            )

    return deduplicate_alternatives(alternatives)


def load_3seq(path: Path) -> Tuple[List[ThreeSeqEvent], List[str]]:
    """Load normalized 3SEQ events or a native ``.3s.rec`` table."""

    rows, headers, _ = read_table(path)
    events: List[ThreeSeqEvent] = []
    warnings: List[str] = []

    for row_number, row in enumerate(rows, start=2):
        recombinant_raw = get_first(
            row,
            "recombinant",
            "recombinant_accession",
            "recombinant_accessions",
            "c_accnum",
            "child",
            "child_accession",
        )
        recombinant = first_accession(recombinant_raw)
        if not recombinant:
            warnings.append(f"3SEQ row {row_number}: missing recombinant identifier; row skipped")
            continue

        parent_1 = first_accession(
            get_first(row, "best_parent_1", "parent_1", "p_accnum", "p_accession")
        )
        parent_2 = first_accession(
            get_first(row, "best_parent_2", "parent_2", "q_accnum", "q_accession")
        )
        event_id = get_first(row, "event_id", "event", "three_seq_event_id")
        if not event_id:
            event_id = f"3SEQ_ROW_{row_number:06d}"

        alternatives = extract_3seq_alternatives(row)
        if not alternatives:
            warnings.append(
                f"3SEQ row {row_number} ({event_id}): no parseable breakpoint pair"
            )

        events.append(
            ThreeSeqEvent(
                event_id=event_id,
                recombinant=recombinant,
                parent_1=parent_1,
                parent_2=parent_2,
                alternatives=alternatives,
                p_value=parse_float(get_first(row, "p_value", "p", "pvalue")),
                min_recombinant_length=parse_int(
                    get_first(row, "min_recombinant_length", "min_rec_length")
                ),
                source_row=row_number,
            )
        )

    if not events:
        raise ValueError(f"No usable 3SEQ events found in {path}")
    return events, warnings


RDP_COORDINATE_COLUMNS = {
    "alignment": ("alignment_begin", "alignment_end", "crosses_alignment_origin"),
    "reference": ("reference_begin", "reference_end", "crosses_reference_origin"),
    "recombinant": ("recombinant_begin", "recombinant_end", ""),
}


def load_rdp(
    path: Path,
    coordinate_system: str,
    include_nonrobust: bool,
) -> Tuple[List[RDPEvent], List[int], List[str]]:
    """Load RDP events and expand rows containing multiple recombinants."""

    rows, _, _ = read_table(path)
    events: List[RDPEvent] = []
    lengths: List[int] = []
    warnings: List[str] = []
    begin_column, end_column, crossing_column = RDP_COORDINATE_COLUMNS[coordinate_system]

    for row_number, row in enumerate(rows, start=2):
        robust_value = get_first(row, "robust_event")
        robust = parse_bool(robust_value)
        if not include_nonrobust and robust is False:
            continue

        recombinant_raw = get_first(
            row,
            "recombinant_accessions",
            "recombinant_accession",
            "recombinants",
            "recombinant_sequences_raw",
            "recombinant_sequences",
            "recombinant_sequence_s",
        )
        recombinants = extract_accessions(recombinant_raw)
        if not recombinants:
            warnings.append(f"RDP row {row_number}: missing recombinant identifier; row skipped")
            continue

        event_id = get_first(row, "event", "event_id", "recombination_event_number")
        if not event_id:
            event_id = f"RDP_ROW_{row_number:06d}"

        minor_parents = tuple(
            extract_accessions(
                get_first(
                    row,
                    "minor_parent_accessions",
                    "minor_parent_candidates_raw",
                    "minor_parent_candidates",
                    "minor_parental_sequence_s",
                )
            )
        )
        major_parents = tuple(
            extract_accessions(
                get_first(
                    row,
                    "major_parent_accessions",
                    "major_parent_candidates_raw",
                    "major_parent_candidates",
                    "major_parental_sequence_s",
                )
            )
        )

        breakpoint_1 = parse_int(get_first(row, begin_column))
        breakpoint_2 = parse_int(get_first(row, end_column))
        crosses_origin = parse_bool(get_first(row, crossing_column)) if crossing_column else None
        p_value = parse_float(get_first(row, "minimum_pvalue", "minimum_p_value"))

        length = parse_int(
            get_first(
                row,
                f"{coordinate_system}_length",
                "alignment_length" if coordinate_system == "alignment" else "",
                "reference_genome_length" if coordinate_system == "reference" else "",
            )
        )
        if length:
            lengths.append(length)

        if breakpoint_1 is None or breakpoint_2 is None:
            warnings.append(
                f"RDP row {row_number} ({event_id}): incomplete {coordinate_system} breakpoints"
            )

        for recombinant in recombinants:
            events.append(
                RDPEvent(
                    event_id=event_id,
                    recombinant=recombinant,
                    minor_parents=minor_parents,
                    major_parents=major_parents,
                    breakpoint_1=breakpoint_1,
                    breakpoint_2=breakpoint_2,
                    crosses_origin=crosses_origin,
                    p_value=p_value,
                    source_row=row_number,
                )
            )

    if not events:
        raise ValueError(f"No usable RDP events found in {path}")
    return events, lengths, warnings


def alignment_length_from_fasta(path: Path) -> int:
    """Return the common sequence length of an aligned FASTA file."""

    lengths: List[int] = []
    current: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(sum(len(part) for part in current))
                    current = []
            else:
                current.append(re.sub(r"\s+", "", line))
        if current:
            lengths.append(sum(len(part) for part in current))

    if not lengths:
        raise ValueError(f"No sequences found in alignment: {path}")
    distinct = sorted(set(lengths))
    if len(distinct) != 1:
        raise ValueError(
            f"FASTA is not aligned: observed sequence lengths are {distinct[:10]}"
        )
    return distinct[0]


def resolve_coordinate_length(
    explicit_length: Optional[int],
    alignment_path: Optional[Path],
    rdp_lengths: Sequence[int],
    topology: str,
) -> Optional[int]:
    """Resolve the coordinate-system length and detect conflicting evidence."""

    fasta_length = alignment_length_from_fasta(alignment_path) if alignment_path else None
    distinct_rdp_lengths = sorted(set(length for length in rdp_lengths if length > 0))
    rdp_length = distinct_rdp_lengths[0] if len(distinct_rdp_lengths) == 1 else None

    if len(distinct_rdp_lengths) > 1:
        eprint(
            "WARNING: several RDP coordinate lengths were found:",
            ", ".join(map(str, distinct_rdp_lengths)),
        )

    if explicit_length is not None:
        if explicit_length <= 0:
            raise ValueError("--coordinate-length must be greater than zero")
        chosen = explicit_length
    elif fasta_length is not None:
        chosen = fasta_length
    else:
        chosen = rdp_length

    if chosen is not None and fasta_length is not None and chosen != fasta_length:
        raise ValueError(
            f"Coordinate length ({chosen}) differs from FASTA alignment length ({fasta_length})"
        )
    if chosen is not None and rdp_length is not None and chosen != rdp_length:
        raise ValueError(
            f"Coordinate length ({chosen}) differs from RDP alignment length ({rdp_length})"
        )
    if topology == "circular" and chosen is None:
        raise ValueError(
            "Circular comparison requires --coordinate-length, --alignment, "
            "or a consistent alignment_length column in the RDP table"
        )
    return chosen


def normalize_circular_position(position: float, length: int) -> float:
    """Map an arbitrary coordinate to the half-open circular range [0, length)."""

    return position % length


def position_to_zero_based(position: int) -> float:
    """Convert a one-based breakpoint coordinate to a zero-based location."""

    return float(position - 1)


def circular_distance(a: float, b: float, length: int) -> float:
    """Shortest distance between two locations on a circle."""

    a = normalize_circular_position(a, length)
    b = normalize_circular_position(b, length)
    direct = abs(a - b)
    return min(direct, length - direct)


def interval_contains_point(
    interval: Interval,
    point: int,
    topology: str,
    length: Optional[int],
) -> bool:
    """Test whether a one-based point lies in an interval."""

    if topology == "linear" or interval.start <= interval.end:
        return interval.start <= point <= interval.end
    if length is None:
        return False
    return point >= interval.start or point <= interval.end


def interval_point_distance(
    interval: Interval,
    point: int,
    topology: str,
    length: Optional[int],
) -> float:
    """Minimum distance between a breakpoint range and an RDP breakpoint."""

    if interval_contains_point(interval, point, topology, length):
        return 0.0

    if topology == "circular":
        if length is None:
            raise ValueError("Circular distance requires a coordinate length")
        point0 = position_to_zero_based(point)
        return min(
            circular_distance(position_to_zero_based(interval.start), point0, length),
            circular_distance(position_to_zero_based(interval.end), point0, length),
        )

    if interval.start > interval.end:
        raise ValueError(
            f"Wrapped interval {interval.label()} cannot be used with linear topology"
        )
    if point < interval.start:
        return float(interval.start - point)
    return float(point - interval.end)


def interval_midpoint(
    interval: Interval,
    topology: str,
    length: Optional[int],
) -> float:
    """Return a representative zero-based midpoint for a breakpoint range."""

    start = position_to_zero_based(interval.start)
    end = position_to_zero_based(interval.end)
    if topology == "circular" and interval.start > interval.end:
        if length is None:
            raise ValueError("Circular midpoint requires a coordinate length")
        span = (end - start) % length
        return normalize_circular_position(start + span / 2.0, length)
    return (start + end) / 2.0


def arc_segments(start: float, end: float, length: float) -> List[Tuple[float, float]]:
    """Represent the clockwise half-open arc [start, end) as linear segments."""

    start %= length
    end %= length
    if math.isclose(start, end):
        return []
    if start < end:
        return [(start, end)]
    return [(start, length), (0.0, end)]


def segments_length(segments: Sequence[Tuple[float, float]]) -> float:
    return sum(max(0.0, end - start) for start, end in segments)


def segments_intersection_length(
    first: Sequence[Tuple[float, float]],
    second: Sequence[Tuple[float, float]],
) -> float:
    total = 0.0
    for start_a, end_a in first:
        for start_b, end_b in second:
            total += max(0.0, min(end_a, end_b) - max(start_a, start_b))
    return total


def segment_jaccard(
    first: Sequence[Tuple[float, float]],
    second: Sequence[Tuple[float, float]],
) -> float:
    intersection = segments_intersection_length(first, second)
    union = segments_length(first) + segments_length(second) - intersection
    return intersection / union if union > 0 else 1.0


def breakpoint_partition_jaccard(
    alternative: BreakpointAlternative,
    rdp_breakpoint_1: int,
    rdp_breakpoint_2: int,
    topology: str,
    length: Optional[int],
) -> Optional[float]:
    """Compare the genome partitions induced by two breakpoint pairs.

    For circular genomes, each breakpoint pair defines two complementary arcs.
    The maximum Jaccard index across both possible arc orientations is returned,
    because standalone 3SEQ does not necessarily label the same recombinant arc
    orientation as RDP.
    """

    if topology != "circular" or length is None:
        return None

    a1 = interval_midpoint(alternative.breakpoint_1, topology, length)
    a2 = interval_midpoint(alternative.breakpoint_2, topology, length)
    b1 = position_to_zero_based(rdp_breakpoint_1) % length
    b2 = position_to_zero_based(rdp_breakpoint_2) % length

    a_forward = arc_segments(a1, a2, float(length))
    a_reverse = arc_segments(a2, a1, float(length))
    b_forward = arc_segments(b1, b2, float(length))
    b_reverse = arc_segments(b2, b1, float(length))

    return max(
        segment_jaccard(a_forward, b_forward),
        segment_jaccard(a_forward, b_reverse),
        segment_jaccard(a_reverse, b_forward),
        segment_jaccard(a_reverse, b_reverse),
    )


def compare_alternatives(
    threeseq: ThreeSeqEvent,
    rdp: RDPEvent,
    topology: str,
    length: Optional[int],
) -> Optional[AlternativeComparison]:
    """Choose the 3SEQ alternative and breakpoint order closest to an RDP event."""

    if (
        not threeseq.alternatives
        or rdp.breakpoint_1 is None
        or rdp.breakpoint_2 is None
    ):
        return None

    candidates: List[AlternativeComparison] = []
    for alternative_index, alternative in enumerate(threeseq.alternatives, start=1):
        for orientation, target_1, target_2 in (
            ("direct", rdp.breakpoint_1, rdp.breakpoint_2),
            ("swapped", rdp.breakpoint_2, rdp.breakpoint_1),
        ):
            distance_1 = interval_point_distance(
                alternative.breakpoint_1, target_1, topology, length
            )
            distance_2 = interval_point_distance(
                alternative.breakpoint_2, target_2, topology, length
            )
            candidates.append(
                AlternativeComparison(
                    alternative_index=alternative_index,
                    alternative=alternative,
                    orientation=orientation,
                    matched_rdp_for_bp1=target_1,
                    matched_rdp_for_bp2=target_2,
                    distance_1=distance_1,
                    distance_2=distance_2,
                    maximum_distance=max(distance_1, distance_2),
                    sum_distance=distance_1 + distance_2,
                    partition_jaccard=breakpoint_partition_jaccard(
                        alternative,
                        rdp.breakpoint_1,
                        rdp.breakpoint_2,
                        topology,
                        length,
                    ),
                )
            )

    return min(
        candidates,
        key=lambda item: (
            item.maximum_distance if item.maximum_distance is not None else math.inf,
            item.sum_distance if item.sum_distance is not None else math.inf,
            item.alternative_index,
            0 if item.orientation == "direct" else 1,
        ),
    )


def assess_parent_concordance(threeseq: ThreeSeqEvent, rdp: RDPEvent) -> str:
    """Assess parent agreement without imposing minor/major orientation."""

    three_parents = {p for p in (threeseq.parent_1, threeseq.parent_2) if p}
    rdp_parents = set(rdp.minor_parents) | set(rdp.major_parents)
    if not rdp_parents:
        return "Unknown_RDP_parent"
    if not three_parents:
        return "Unknown_3SEQ_parent"

    overlap = three_parents & rdp_parents
    if len(three_parents) >= 2 and three_parents.issubset(rdp_parents):
        return "Both_parents"
    if overlap:
        return "One_parent"
    return "No_parent"


def classify_concordance(
    best: Optional[AlternativeComparison],
    tolerance: float,
    moderate_multiplier: float,
) -> str:
    """Assign a breakpoint concordance class."""

    if best is None or best.distance_1 is None or best.distance_2 is None:
        return "Insufficient_breakpoint_data"
    if best.distance_1 <= tolerance and best.distance_2 <= tolerance:
        return "Strong_concordance"
    if (
        best.distance_1 <= tolerance * moderate_multiplier
        and best.distance_2 <= tolerance * moderate_multiplier
    ):
        return "Moderate_concordance"
    if min(best.distance_1, best.distance_2) <= tolerance:
        return "Partial_concordance"
    return "Sequence_only"


def event_pair_cost(
    threeseq: ThreeSeqEvent,
    rdp: RDPEvent,
    topology: str,
    length: Optional[int],
    parent_penalty: float,
) -> Tuple[float, Optional[AlternativeComparison], str]:
    """Calculate assignment cost and retain the best breakpoint comparison."""

    best = compare_alternatives(threeseq, rdp, topology, length)
    parent_concordance = assess_parent_concordance(threeseq, rdp)

    if best is None or best.maximum_distance is None or best.sum_distance is None:
        base_cost = 1.0e12
    else:
        # Maximum distance is primary; the small sum term breaks deterministic ties.
        base_cost = best.maximum_distance + best.sum_distance * 1.0e-6

    if parent_penalty > 0:
        mismatch_units = {
            "Both_parents": 0,
            "One_parent": 1,
            "No_parent": 2,
            "Unknown_RDP_parent": 0,
            "Unknown_3SEQ_parent": 0,
        }[parent_concordance]
        base_cost += parent_penalty * mismatch_units

    return base_cost, best, parent_concordance


def linear_sum_assignment(cost_matrix: Sequence[Sequence[float]]) -> List[Tuple[int, int]]:
    """Solve a rectangular minimum-cost assignment using the Hungarian method.

    This standard-library implementation returns ``min(rows, columns)`` pairs.
    It avoids requiring SciPy in the Snakemake environment.
    """

    if not cost_matrix or not cost_matrix[0]:
        return []
    rows = len(cost_matrix)
    columns = len(cost_matrix[0])
    if any(len(row) != columns for row in cost_matrix):
        raise ValueError("Cost matrix is not rectangular")

    transposed = rows > columns
    matrix = (
        [[cost_matrix[i][j] for i in range(rows)] for j in range(columns)]
        if transposed
        else [list(row) for row in cost_matrix]
    )
    n = len(matrix)
    m = len(matrix[0])

    u = [0.0] * (n + 1)
    v = [0.0] * (m + 1)
    p = [0] * (m + 1)
    way = [0] * (m + 1)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minimum = [math.inf] * (m + 1)
        used = [False] * (m + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = math.inf
            j1 = 0
            for j in range(1, m + 1):
                if used[j]:
                    continue
                current = matrix[i0 - 1][j - 1] - u[i0] - v[j]
                if current < minimum[j]:
                    minimum[j] = current
                    way[j] = j0
                if minimum[j] < delta:
                    delta = minimum[j]
                    j1 = j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minimum[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    assignments: List[Tuple[int, int]] = []
    for column_index in range(1, m + 1):
        if p[column_index] == 0:
            continue
        row_index = p[column_index] - 1
        col_index = column_index - 1
        assignments.append(
            (col_index, row_index) if transposed else (row_index, col_index)
        )
    return sorted(assignments)


def compare_events(
    threeseq_events: Sequence[ThreeSeqEvent],
    rdp_events: Sequence[RDPEvent],
    topology: str,
    length: Optional[int],
    tolerance: float,
    moderate_multiplier: float,
    parent_penalty: float,
) -> List[EventComparison]:
    """Perform accession-stratified, one-to-one event matching."""

    three_by_recombinant: Dict[str, List[ThreeSeqEvent]] = defaultdict(list)
    rdp_by_recombinant: Dict[str, List[RDPEvent]] = defaultdict(list)
    for event in threeseq_events:
        three_by_recombinant[event.recombinant].append(event)
    for event in rdp_events:
        rdp_by_recombinant[event.recombinant].append(event)

    comparisons: List[EventComparison] = []
    for recombinant in sorted(set(three_by_recombinant) | set(rdp_by_recombinant)):
        three_group = sorted(
            three_by_recombinant.get(recombinant, []),
            key=lambda event: (event.event_id, event.source_row),
        )
        rdp_group = sorted(
            rdp_by_recombinant.get(recombinant, []),
            key=lambda event: (event.event_id, event.source_row),
        )

        if not three_group:
            for rdp in rdp_group:
                comparisons.append(
                    EventComparison(
                        recombinant,
                        None,
                        rdp,
                        "RDP_only",
                        "RDP_only",
                        "Not_applicable",
                        None,
                        None,
                    )
                )
            continue
        if not rdp_group:
            for threeseq in three_group:
                comparisons.append(
                    EventComparison(
                        recombinant,
                        threeseq,
                        None,
                        "3SEQ_only",
                        "3SEQ_only",
                        "Not_applicable",
                        None,
                        None,
                    )
                )
            continue

        cost_matrix: List[List[float]] = []
        pair_details: Dict[Tuple[int, int], Tuple[Optional[AlternativeComparison], str]] = {}
        for i, threeseq in enumerate(three_group):
            row_costs: List[float] = []
            for j, rdp in enumerate(rdp_group):
                cost, best, parent_concordance = event_pair_cost(
                    threeseq,
                    rdp,
                    topology,
                    length,
                    parent_penalty,
                )
                row_costs.append(cost)
                pair_details[(i, j)] = (best, parent_concordance)
            cost_matrix.append(row_costs)

        assignments = linear_sum_assignment(cost_matrix)
        matched_three = set()
        matched_rdp = set()
        for i, j in assignments:
            matched_three.add(i)
            matched_rdp.add(j)
            best, parent_concordance = pair_details[(i, j)]
            comparisons.append(
                EventComparison(
                    recombinant=recombinant,
                    threeseq=three_group[i],
                    rdp=rdp_group[j],
                    match_status="Matched",
                    concordance_class=classify_concordance(
                        best, tolerance, moderate_multiplier
                    ),
                    parent_concordance=parent_concordance,
                    best=best,
                    assignment_cost=cost_matrix[i][j],
                )
            )

        for i, threeseq in enumerate(three_group):
            if i not in matched_three:
                comparisons.append(
                    EventComparison(
                        recombinant,
                        threeseq,
                        None,
                        "3SEQ_only",
                        "3SEQ_only",
                        "Not_applicable",
                        None,
                        None,
                    )
                )
        for j, rdp in enumerate(rdp_group):
            if j not in matched_rdp:
                comparisons.append(
                    EventComparison(
                        recombinant,
                        None,
                        rdp,
                        "RDP_only",
                        "RDP_only",
                        "Not_applicable",
                        None,
                        None,
                    )
                )

    return sorted(
        comparisons,
        key=lambda item: (
            item.recombinant,
            item.threeseq.event_id if item.threeseq else "~",
            item.rdp.event_id if item.rdp else "~",
        ),
    )


PAIR_FIELDS = [
    "Recombinant_accession",
    "Match_status",
    "Concordance_class",
    "3SEQ_Event_ID",
    "RDP_Event_ID",
    "3SEQ_source_row",
    "RDP_source_row",
    "3SEQ_parent_1",
    "3SEQ_parent_2",
    "RDP_minor_parents",
    "RDP_major_parents",
    "Parent_concordance",
    "3SEQ_breakpoint_alternatives_count",
    "Selected_3SEQ_alternative_index",
    "Selected_3SEQ_alternative_source",
    "Selected_3SEQ_breakpoint_pair",
    "3SEQ_BP1_start",
    "3SEQ_BP1_end",
    "3SEQ_BP2_start",
    "3SEQ_BP2_end",
    "RDP_breakpoint_1",
    "RDP_breakpoint_2",
    "Breakpoint_order",
    "RDP_breakpoint_matched_to_3SEQ_BP1",
    "RDP_breakpoint_matched_to_3SEQ_BP2",
    "Breakpoint_1_distance_nt",
    "Breakpoint_2_distance_nt",
    "Maximum_breakpoint_distance_nt",
    "Sum_breakpoint_distance_nt",
    "Breakpoint_partition_Jaccard",
    "3SEQ_P_value",
    "RDP_minimum_P_value",
    "3SEQ_min_recombinant_length",
    "RDP_crosses_origin",
    "Assignment_cost",
]


def comparison_to_row(comparison: EventComparison) -> Dict[str, object]:
    """Flatten an event comparison for tabular output."""

    threeseq = comparison.threeseq
    rdp = comparison.rdp
    best = comparison.best
    alternative = best.alternative if best else None

    return {
        "Recombinant_accession": comparison.recombinant,
        "Match_status": comparison.match_status,
        "Concordance_class": comparison.concordance_class,
        "3SEQ_Event_ID": threeseq.event_id if threeseq else "",
        "RDP_Event_ID": rdp.event_id if rdp else "",
        "3SEQ_source_row": threeseq.source_row if threeseq else "",
        "RDP_source_row": rdp.source_row if rdp else "",
        "3SEQ_parent_1": threeseq.parent_1 if threeseq else "",
        "3SEQ_parent_2": threeseq.parent_2 if threeseq else "",
        "RDP_minor_parents": ";".join(rdp.minor_parents) if rdp else "",
        "RDP_major_parents": ";".join(rdp.major_parents) if rdp else "",
        "Parent_concordance": comparison.parent_concordance,
        "3SEQ_breakpoint_alternatives_count": len(threeseq.alternatives) if threeseq else "",
        "Selected_3SEQ_alternative_index": best.alternative_index if best else "",
        "Selected_3SEQ_alternative_source": alternative.source if alternative else "",
        "Selected_3SEQ_breakpoint_pair": alternative.label() if alternative else "",
        "3SEQ_BP1_start": alternative.breakpoint_1.start if alternative else "",
        "3SEQ_BP1_end": alternative.breakpoint_1.end if alternative else "",
        "3SEQ_BP2_start": alternative.breakpoint_2.start if alternative else "",
        "3SEQ_BP2_end": alternative.breakpoint_2.end if alternative else "",
        "RDP_breakpoint_1": rdp.breakpoint_1 if rdp and rdp.breakpoint_1 is not None else "",
        "RDP_breakpoint_2": rdp.breakpoint_2 if rdp and rdp.breakpoint_2 is not None else "",
        "Breakpoint_order": best.orientation if best else "",
        "RDP_breakpoint_matched_to_3SEQ_BP1": best.matched_rdp_for_bp1 if best else "",
        "RDP_breakpoint_matched_to_3SEQ_BP2": best.matched_rdp_for_bp2 if best else "",
        "Breakpoint_1_distance_nt": format_number(best.distance_1, 3) if best else "",
        "Breakpoint_2_distance_nt": format_number(best.distance_2, 3) if best else "",
        "Maximum_breakpoint_distance_nt": (
            format_number(best.maximum_distance, 3) if best else ""
        ),
        "Sum_breakpoint_distance_nt": format_number(best.sum_distance, 3) if best else "",
        "Breakpoint_partition_Jaccard": (
            format_number(best.partition_jaccard, 6) if best else ""
        ),
        "3SEQ_P_value": (
            format_number(threeseq.p_value, 12) if threeseq and threeseq.p_value is not None else ""
        ),
        "RDP_minimum_P_value": (
            format_number(rdp.p_value, 12) if rdp and rdp.p_value is not None else ""
        ),
        "3SEQ_min_recombinant_length": (
            threeseq.min_recombinant_length
            if threeseq and threeseq.min_recombinant_length is not None
            else ""
        ),
        "RDP_crosses_origin": (
            str(rdp.crosses_origin) if rdp and rdp.crosses_origin is not None else ""
        ),
        "Assignment_cost": format_number(comparison.assignment_cost, 6),
    }


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    """Write a deterministic UTF-8 TSV file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_recombinant_summary(
    comparisons: Sequence[EventComparison],
) -> List[Dict[str, object]]:
    """Summarize matched and unmatched events for each recombinant accession."""

    by_recombinant: Dict[str, List[EventComparison]] = defaultdict(list)
    for comparison in comparisons:
        by_recombinant[comparison.recombinant].append(comparison)

    result: List[Dict[str, object]] = []
    for recombinant in sorted(by_recombinant):
        group = by_recombinant[recombinant]
        three_ids = {c.threeseq.event_id for c in group if c.threeseq}
        rdp_ids = {c.rdp.event_id for c in group if c.rdp}
        classes = Counter(c.concordance_class for c in group)
        parents = Counter(c.parent_concordance for c in group if c.match_status == "Matched")
        distances = [
            c.best.maximum_distance
            for c in group
            if c.best is not None and c.best.maximum_distance is not None
        ]
        jaccards = [
            c.best.partition_jaccard
            for c in group
            if c.best is not None and c.best.partition_jaccard is not None
        ]

        result.append(
            {
                "Recombinant_accession": recombinant,
                "N_3SEQ_events": len(three_ids),
                "N_RDP_events": len(rdp_ids),
                "N_matched_event_pairs": sum(c.match_status == "Matched" for c in group),
                "N_strong_concordance": classes["Strong_concordance"],
                "N_moderate_concordance": classes["Moderate_concordance"],
                "N_partial_concordance": classes["Partial_concordance"],
                "N_sequence_only": classes["Sequence_only"],
                "N_insufficient_breakpoint_data": classes["Insufficient_breakpoint_data"],
                "N_3SEQ_only_events": classes["3SEQ_only"],
                "N_RDP_only_events": classes["RDP_only"],
                "Best_maximum_breakpoint_distance_nt": (
                    format_number(min(distances), 3) if distances else ""
                ),
                "Best_breakpoint_partition_Jaccard": (
                    format_number(max(jaccards), 6) if jaccards else ""
                ),
                "N_both_parents": parents["Both_parents"],
                "N_one_parent": parents["One_parent"],
                "N_no_parent": parents["No_parent"],
            }
        )
    return result


RECOMBINANT_SUMMARY_FIELDS = [
    "Recombinant_accession",
    "N_3SEQ_events",
    "N_RDP_events",
    "N_matched_event_pairs",
    "N_strong_concordance",
    "N_moderate_concordance",
    "N_partial_concordance",
    "N_sequence_only",
    "N_insufficient_breakpoint_data",
    "N_3SEQ_only_events",
    "N_RDP_only_events",
    "Best_maximum_breakpoint_distance_nt",
    "Best_breakpoint_partition_Jaccard",
    "N_both_parents",
    "N_one_parent",
    "N_no_parent",
]


def build_global_summary(
    comparisons: Sequence[EventComparison],
    threeseq_events: Sequence[ThreeSeqEvent],
    rdp_events: Sequence[RDPEvent],
    topology: str,
    length: Optional[int],
    tolerance: float,
    moderate_multiplier: float,
    coordinate_system: str,
) -> List[Dict[str, object]]:
    """Build global input and concordance metrics."""

    class_counts = Counter(c.concordance_class for c in comparisons)
    parent_counts = Counter(
        c.parent_concordance for c in comparisons if c.match_status == "Matched"
    )
    matched = sum(c.match_status == "Matched" for c in comparisons)
    strong_accessions = {
        c.recombinant for c in comparisons if c.concordance_class == "Strong_concordance"
    }
    threeseq_accessions = {e.recombinant for e in threeseq_events}
    rdp_accessions = {e.recombinant for e in rdp_events}
    shared_accessions = threeseq_accessions & rdp_accessions

    metrics: List[Tuple[str, object]] = [
        ("Topology", topology),
        ("RDP_coordinate_system", coordinate_system),
        ("Coordinate_length", length if length is not None else ""),
        ("Strong_tolerance_nt", tolerance),
        ("Moderate_tolerance_multiplier", moderate_multiplier),
        ("Input_3SEQ_events", len(threeseq_events)),
        ("Input_RDP_events_after_recombinant_expansion", len(rdp_events)),
        ("Unique_3SEQ_recombinants", len(threeseq_accessions)),
        ("Unique_RDP_recombinants", len(rdp_accessions)),
        ("Shared_recombinant_accessions", len(shared_accessions)),
        ("Matched_event_pairs", matched),
        ("Recombinants_with_strong_event_concordance", len(strong_accessions)),
    ]
    for class_name in (
        "Strong_concordance",
        "Moderate_concordance",
        "Partial_concordance",
        "Sequence_only",
        "Insufficient_breakpoint_data",
        "3SEQ_only",
        "RDP_only",
    ):
        metrics.append((f"Event_class_{class_name}", class_counts[class_name]))
    for parent_class in (
        "Both_parents",
        "One_parent",
        "No_parent",
        "Unknown_RDP_parent",
        "Unknown_3SEQ_parent",
    ):
        metrics.append((f"Parent_class_{parent_class}", parent_counts[parent_class]))

    if matched:
        metrics.append(
            (
                "Strong_event_concordance_percent_of_matched_pairs",
                f"{100.0 * class_counts['Strong_concordance'] / matched:.3f}",
            )
        )
    if rdp_accessions:
        metrics.append(
            (
                "RDP_recombinants_detected_by_3SEQ_percent",
                f"{100.0 * len(shared_accessions) / len(rdp_accessions):.3f}",
            )
        )
    if shared_accessions:
        metrics.append(
            (
                "Shared_accessions_with_strong_breakpoint_concordance_percent",
                f"{100.0 * len(strong_accessions) / len(shared_accessions):.3f}",
            )
        )

    return [{"Metric": metric, "Value": value} for metric, value in metrics]


def validate_event_coordinates(
    threeseq_events: Sequence[ThreeSeqEvent],
    rdp_events: Sequence[RDPEvent],
    length: Optional[int],
    topology: str,
) -> List[str]:
    """Report coordinates outside the declared coordinate system.

    Standalone 3SEQ can report ``0`` as the breakpoint boundary immediately
    before the first alignment site.  This is a valid origin boundary for a
    circular alignment and is therefore accepted only for circular 3SEQ data.
    RDP alignment coordinates remain one-based.
    """

    if length is None:
        return []
    warnings: List[str] = []
    observed: List[Tuple[str, str, int]] = []
    for event in threeseq_events:
        for alternative in event.alternatives:
            observed.extend(
                [
                    ("3SEQ", event.event_id, alternative.breakpoint_1.start),
                    ("3SEQ", event.event_id, alternative.breakpoint_1.end),
                    ("3SEQ", event.event_id, alternative.breakpoint_2.start),
                    ("3SEQ", event.event_id, alternative.breakpoint_2.end),
                ]
            )
    for event in rdp_events:
        if event.breakpoint_1 is not None:
            observed.append(("RDP", event.event_id, event.breakpoint_1))
        if event.breakpoint_2 is not None:
            observed.append(("RDP", event.event_id, event.breakpoint_2))

    for source, event_id, coordinate in observed:
        minimum = 0 if source == "3SEQ" and topology == "circular" else 1
        if coordinate < minimum or coordinate > length:
            warnings.append(
                f"{source} event {event_id}: coordinate {coordinate} is outside "
                f"{minimum}..{length}"
            )
    return warnings


def output_paths(prefix: Path) -> Tuple[Path, Path, Path]:
    """Derive output files from a user-supplied prefix."""

    text = str(prefix)
    if text.endswith(".tsv"):
        text = text[:-4]
    base = Path(text)
    return (
        Path(f"{base}.event_pairs.tsv"),
        Path(f"{base}.recombinant_summary.tsv"),
        Path(f"{base}.summary.tsv"),
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare standalone 3SEQ and RDP recombination events using "
            "one-to-one breakpoint matching."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--threeseq", required=True, type=Path, help="3SEQ TSV or .3s.rec")
    parser.add_argument("--rdp", required=True, type=Path, help="Parsed RDP event TSV")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output file prefix")
    parser.add_argument(
        "--alignment",
        type=Path,
        help="Aligned FASTA used by both programs; supplies the coordinate length",
    )
    parser.add_argument(
        "--coordinate-length",
        "--genome-length",
        dest="coordinate_length",
        type=int,
        help="Length of the selected coordinate system",
    )
    parser.add_argument(
        "--topology",
        choices=("circular", "linear"),
        default="circular",
        help="Genome/alignment topology used for breakpoint distances",
    )
    parser.add_argument(
        "--rdp-coordinate-system",
        choices=tuple(RDP_COORDINATE_COLUMNS),
        default="alignment",
        help="RDP coordinate columns to compare with 3SEQ",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=50.0,
        help="Maximum distance for strong concordance of each breakpoint",
    )
    parser.add_argument(
        "--moderate-multiplier",
        type=float,
        default=2.0,
        help="Moderate threshold as a multiple of --tolerance",
    )
    parser.add_argument(
        "--parent-penalty",
        type=float,
        default=0.0,
        help=(
            "Optional cost added for parent mismatch during event assignment; "
            "zero keeps breakpoint concordance primary"
        ),
    )
    parser.add_argument(
        "--include-nonrobust-rdp",
        action="store_true",
        help="Retain RDP rows explicitly marked robust_event=False",
    )
    parser.add_argument(
        "--max-warnings",
        type=int,
        default=30,
        help="Maximum number of individual parsing warnings printed",
    )
    args = parser.parse_args(argv)

    if args.tolerance < 0:
        parser.error("--tolerance must be non-negative")
    if args.moderate_multiplier < 1:
        parser.error("--moderate-multiplier must be at least 1")
    if args.parent_penalty < 0:
        parser.error("--parent-penalty must be non-negative")
    if args.max_warnings < 0:
        parser.error("--max-warnings must be non-negative")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    try:
        if not args.threeseq.is_file():
            raise FileNotFoundError(f"3SEQ input not found: {args.threeseq}")
        if not args.rdp.is_file():
            raise FileNotFoundError(f"RDP input not found: {args.rdp}")
        if args.alignment is not None and not args.alignment.is_file():
            raise FileNotFoundError(f"Alignment not found: {args.alignment}")

        threeseq_events, three_warnings = load_3seq(args.threeseq)
        rdp_events, rdp_lengths, rdp_warnings = load_rdp(
            args.rdp,
            args.rdp_coordinate_system,
            args.include_nonrobust_rdp,
        )
        length = resolve_coordinate_length(
            args.coordinate_length,
            args.alignment,
            rdp_lengths,
            args.topology,
        )
        coordinate_warnings = validate_event_coordinates(
            threeseq_events, rdp_events, length, args.topology
        )

        warnings = three_warnings + rdp_warnings + coordinate_warnings
        for warning in warnings[: args.max_warnings]:
            eprint("WARNING:", warning)
        if len(warnings) > args.max_warnings:
            eprint(
                f"WARNING: {len(warnings) - args.max_warnings} additional warnings suppressed"
            )

        comparisons = compare_events(
            threeseq_events,
            rdp_events,
            args.topology,
            length,
            args.tolerance,
            args.moderate_multiplier,
            args.parent_penalty,
        )

        pair_path, recombinant_path, summary_path = output_paths(args.out_prefix)
        write_tsv(pair_path, PAIR_FIELDS, map(comparison_to_row, comparisons))
        write_tsv(
            recombinant_path,
            RECOMBINANT_SUMMARY_FIELDS,
            build_recombinant_summary(comparisons),
        )
        write_tsv(
            summary_path,
            ["Metric", "Value"],
            build_global_summary(
                comparisons,
                threeseq_events,
                rdp_events,
                args.topology,
                length,
                args.tolerance,
                args.moderate_multiplier,
                args.rdp_coordinate_system,
            ),
        )

        classes = Counter(c.concordance_class for c in comparisons)
        matched = sum(c.match_status == "Matched" for c in comparisons)
        print(f"3SEQ events                         : {len(threeseq_events)}")
        print(f"RDP events (expanded by accession) : {len(rdp_events)}")
        print(f"Matched event pairs                 : {matched}")
        print(f"Strong concordance                  : {classes['Strong_concordance']}")
        print(f"Moderate concordance                : {classes['Moderate_concordance']}")
        print(f"Partial concordance                 : {classes['Partial_concordance']}")
        print(f"Sequence-only matches               : {classes['Sequence_only']}")
        print(f"3SEQ-only events                    : {classes['3SEQ_only']}")
        print(f"RDP-only events                     : {classes['RDP_only']}")
        print(f"Detailed pairs written to           : {pair_path}")
        print(f"Recombinant summary written to      : {recombinant_path}")
        print(f"Global summary written to           : {summary_path}")
        return 0

    except (OSError, ValueError, csv.Error) as error:
        eprint(f"ERROR: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
