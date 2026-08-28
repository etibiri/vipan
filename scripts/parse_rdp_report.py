#!/usr/bin/env python3

import argparse
import csv
import json
import math
import re
from pathlib import Path


CANONICAL_COLUMNS = [
    "event",
    "rdp_file_number",
    "alignment_begin",
    "alignment_end",
    "recombinant_begin",
    "recombinant_end",
    "reference_begin",
    "reference_end",
    "recombinant_sequences",
    "minor_parent_candidates",
    "major_parent_candidates",
    "rdp",
    "geneconv",
    "bootscan",
    "maxchi",
    "chimaera",
    "siscan",
    "phylpro",
    "lard",
    "threeseq",
]

METHOD_COLUMNS = {
    "rdp": 11,
    "geneconv": 12,
    "bootscan": 13,
    "maxchi": 14,
    "chimaera": 15,
    "siscan": 16,
    "phylpro": 17,
    "lard": 18,
    "threeseq": 19,
}

# Méthodes utilisées pour la comparaison indépendante avec 3SEQ.
INDEPENDENT_METHODS = [
    "rdp",
    "geneconv",
    "bootscan",
    "maxchi",
    "chimaera",
    "siscan",
    "phylpro",
    "lard",
]

ALL_METHODS = INDEPENDENT_METHODS + ["threeseq"]

INTEGER_PATTERN = re.compile(r"-?\d+")

PVALUE_PATTERN = re.compile(
    r"([0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)"
)

REFERENCE_PATTERN = re.compile(
    r"Relative\s+to\s+([A-Za-z0-9_.-]+)",
    flags=re.IGNORECASE,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Parse and filter a hierarchical RDP recombination report."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="RDP CSV report.",
    )

    parser.add_argument(
        "--alignment",
        required=True,
        help="Aligned FASTA file used for the RDP analysis.",
    )

    parser.add_argument(
        "--events-all",
        required=True,
        help="Output TSV containing all normalized events.",
    )

    parser.add_argument(
        "--events-robust",
        required=True,
        help="Output TSV containing robust events only.",
    )

    parser.add_argument(
        "--candidates",
        required=True,
        help="Output TSV containing candidate sequence combinations.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Output JSON summary.",
    )

    parser.add_argument(
        "--reference-accession",
        default=None,
        help=(
            "Reference accession. If omitted, it is extracted "
            "from the second header line."
        ),
    )

    parser.add_argument(
        "--min-methods",
        type=int,
        default=3,
        help="Minimum number of significant independent methods.",
    )

    parser.add_argument(
        "--min-tested-methods",
        type=int,
        default=6,
        help="Minimum number of tested independent methods.",
    )

    parser.add_argument(
        "--min-support-fraction",
        type=float,
        default=0.50,
        help="Minimum fraction of significant tested methods.",
    )

    parser.add_argument(
        "--max-missing-methods",
        type=int,
        default=2,
        help="Maximum number of unreported independent methods.",
    )

    parser.add_argument(
        "--method-pvalue-cutoff",
        type=float,
        default=0.05,
        help="Maximum p-value for a method to be significant.",
    )

    parser.add_argument(
        "--min-fragment-length",
        type=int,
        default=100,
        help="Minimum reference recombinant fragment length.",
    )

    parser.add_argument(
        "--allow-caret",
        action="store_true",
        help="Allow events with uncertain recombinant marker (^).",
    )

    parser.add_argument(
        "--allow-tilde",
        action="store_true",
        help="Allow events marked as possibly non-recombinant (~).",
    )

    parser.add_argument(
        "--require-known-parents",
        action="store_true",
        help=(
            "Require at least one known minor and one known "
            "major parent candidate."
        ),
    )

    return parser.parse_args()


def parse_integer(value):
    """Extract an integer while tolerating RDP markers."""
    if value is None:
        return None

    match = INTEGER_PATTERN.search(str(value))

    if match is None:
        return None

    return int(match.group())


def parse_pvalue(value):
    """Extract a numeric p-value from an RDP method cell."""
    if value is None:
        return None

    text = str(value).strip()

    if not text or text.upper() == "NS":
        return None

    match = PVALUE_PATTERN.search(text)

    if match is None:
        return None

    try:
        return float(match.group(1))
    except ValueError:
        return None


def method_status(value, cutoff):
    """
    Classify an RDP method cell.

    Returns
    -------
    status : str
        significant, nonsignificant or not_reported.
    pvalue : float or None
    """
    text = str(value or "").strip()

    if not text:
        return "not_reported", None

    if text.upper() == "NS":
        return "nonsignificant", None

    pvalue = parse_pvalue(text)

    if pvalue is None:
        return "not_reported", None

    if pvalue <= cutoff:
        return "significant", pvalue

    return "nonsignificant", pvalue


def clean_identifier(value):
    """Remove RDP markers while retaining the accession."""
    value = str(value or "").strip()

    if not value:
        return ""

    value = value.replace("^", "")
    value = value.replace("~", "")
    value = value.replace("*", "")
    value = value.strip()

    if value.lower().startswith("unknown"):
        return "Unknown"

    # The accession is expected to be the first token.
    return value.split()[0].strip(";,|")


def unique_preserving_order(values):
    observed = set()
    result = []

    for value in values:
        if not value or value in observed:
            continue

        observed.add(value)
        result.append(value)

    return result


def read_alignment(fasta_file):
    """Read an aligned FASTA file without requiring Biopython."""
    sequences = {}
    current_id = None
    sequence_parts = []

    def save_current():
        if current_id is not None:
            sequences[current_id] = "".join(sequence_parts).upper()

    with open(
        fasta_file,
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_current()

                current_id = line[1:].split()[0]
                sequence_parts = []
            else:
                if current_id is None:
                    raise ValueError(
                        "Sequence encountered before the first FASTA header."
                    )

                sequence_parts.append(line)

    save_current()

    if not sequences:
        raise ValueError(
            f"No sequence found in alignment: {fasta_file}"
        )

    lengths = {len(sequence) for sequence in sequences.values()}

    if len(lengths) != 1:
        raise ValueError(
            "The input FASTA file is not a rectangular alignment."
        )

    alignment_length = next(iter(lengths))

    return sequences, alignment_length


def reference_length_from_alignment(sequences, reference_accession):
    if reference_accession not in sequences:
        raise ValueError(
            f"Reference accession {reference_accession!r} "
            "was not found in the alignment."
        )

    sequence = sequences[reference_accession]

    # N and other nucleotide ambiguity codes still occupy genome positions.
    return sum(
        character not in {"-", "."}
        for character in sequence
    )


def circular_fragment_length(start, end, genome_length):
    if (
        start is None
        or end is None
        or genome_length is None
        or genome_length < 1
    ):
        return None

    if start <= end:
        return end - start + 1

    return genome_length - start + 1 + end


def split_circular_segment(start, end, genome_length):
    """
    Represent an origin-crossing segment as at most two intervals.
    """
    if start is None or end is None:
        return None, None, None, None

    if start <= end:
        return start, end, None, None

    return start, genome_length, 1, end


def empty_event():
    return {
        "event": None,
        "rdp_file_number": None,
        "alignment_begin_raw": "",
        "alignment_end_raw": "",
        "recombinant_begin_raw": "",
        "recombinant_end_raw": "",
        "reference_begin_raw": "",
        "reference_end_raw": "",
        "method_raw": {},
        "recombinant_raw": [],
        "minor_parent_raw": [],
        "major_parent_raw": [],
        "recombinant_accessions": [],
        "minor_parent_accessions": [],
        "major_parent_accessions": [],
        "candidate_rows": [],
        "contains_alternative_process_marker": False,
        "contains_undetermined_breakpoint_marker": False,
        "contains_recombinant_caret_marker": False,
        "contains_parental_marker": False,
        "contains_unknown_minor_parent": False,
        "contains_unknown_major_parent": False,
        "source_line": None,
    }


def add_candidate_row(event, row, line_number):
    recombinant_raw = row[8].strip()
    minor_raw = row[9].strip()
    major_raw = row[10].strip()

    recombinant = clean_identifier(recombinant_raw)
    minor_parent = clean_identifier(minor_raw)
    major_parent = clean_identifier(major_raw)

    if recombinant_raw:
        event["recombinant_raw"].append(recombinant_raw)

    if minor_raw:
        event["minor_parent_raw"].append(minor_raw)

    if major_raw:
        event["major_parent_raw"].append(major_raw)

    if recombinant:
        event["recombinant_accessions"].append(recombinant)

    if minor_parent:
        event["minor_parent_accessions"].append(minor_parent)

    if major_parent:
        event["major_parent_accessions"].append(major_parent)

    candidate_text = " ".join(
        [
            recombinant_raw,
            minor_raw,
            major_raw,
        ]
    )

    event["contains_alternative_process_marker"] |= (
        "~" in candidate_text
    )

    event["contains_recombinant_caret_marker"] |= (
        "^" in recombinant_raw
    )

    event["contains_parental_marker"] |= any(
        marker in (minor_raw + major_raw)
        for marker in ("^", "~", "*")
    )

    event["contains_unknown_minor_parent"] |= (
        minor_parent == "Unknown"
    )

    event["contains_unknown_major_parent"] |= (
        major_parent == "Unknown"
    )

    if recombinant_raw or minor_raw or major_raw:
        event["candidate_rows"].append(
            {
                "source_line": line_number,
                "recombinant_raw": recombinant_raw,
                "recombinant_accession": recombinant,
                "minor_parent_raw": minor_raw,
                "minor_parent_accession": minor_parent,
                "major_parent_raw": major_raw,
                "major_parent_accession": major_parent,
            }
        )


def normalize_row(row, line_number):
    """Normalize main RDP rows and continuation rows to 20 columns."""
    while len(row) > len(CANONICAL_COLUMNS) and not row[-1].strip():
        row.pop()

    if len(row) > len(CANONICAL_COLUMNS):
        raise ValueError(
            f"Line {line_number}: {len(row)} columns after removing "
            "trailing empty cells; expected at most 20."
        )

    if len(row) < len(CANONICAL_COLUMNS):
        row.extend(
            [""] * (len(CANONICAL_COLUMNS) - len(row))
        )

    return row


def parse_rdp_report(report_file, reference_override=None):
    events = []
    current_event = None
    reference_accession = reference_override

    parse_statistics = {
        "blank_lines": 0,
        "metadata_lines": 0,
        "main_event_lines": 0,
        "continuation_lines": 0,
    }

    with open(
        report_file,
        newline="",
        encoding="utf-8-sig",
        errors="replace",
    ) as handle:
        reader = csv.reader(handle)

        header_1 = next(reader, [])
        header_2 = next(reader, [])
        header_3 = next(reader, [])

        if reference_accession is None:
            header_text = " ".join(header_2)
            match = REFERENCE_PATTERN.search(header_text)

            if match is None:
                raise ValueError(
                    "Unable to determine the reference accession "
                    "from the second header line. Use "
                    "--reference-accession."
                )

            reference_accession = match.group(1)

        if len(header_3) != 20:
            raise ValueError(
                "The canonical RDP header does not contain 20 columns."
            )

        for line_number, row in enumerate(reader, start=4):
            if not row or not any(cell.strip() for cell in row):
                parse_statistics["blank_lines"] += 1
                continue

            if (
                len(row) == 1
                and row[0].strip() == "Table Key:"
            ):
                break

            if len(row) == 1:
                parse_statistics["metadata_lines"] += 1
                continue

            row = normalize_row(row, line_number)

            is_main_event = bool(
                row[1].strip()
                and row[2].strip()
                and row[3].strip()
            )

            if is_main_event:
                if current_event is not None:
                    events.append(current_event)

                current_event = empty_event()

                current_event["event"] = parse_integer(row[0])
                current_event["rdp_file_number"] = parse_integer(row[1])
                current_event["source_line"] = line_number

                current_event["alignment_begin_raw"] = row[2].strip()
                current_event["alignment_end_raw"] = row[3].strip()

                current_event["recombinant_begin_raw"] = row[4].strip()
                current_event["recombinant_end_raw"] = row[5].strip()

                current_event["reference_begin_raw"] = row[6].strip()
                current_event["reference_end_raw"] = row[7].strip()

                coordinate_text = " ".join(row[2:8])

                current_event[
                    "contains_undetermined_breakpoint_marker"
                ] = "*" in coordinate_text

                current_event[
                    "contains_alternative_process_marker"
                ] |= "~" in " ".join(row[:11])

                current_event["method_raw"] = {
                    method: row[column].strip()
                    for method, column in METHOD_COLUMNS.items()
                }

                add_candidate_row(
                    current_event,
                    row,
                    line_number,
                )

                parse_statistics["main_event_lines"] += 1

            else:
                if current_event is None:
                    raise ValueError(
                        f"Continuation line {line_number} encountered "
                        "before the first main event."
                    )

                continuation_event = parse_integer(row[0])

                if (
                    continuation_event is not None
                    and continuation_event
                    != current_event["event"]
                ):
                    raise ValueError(
                        f"Line {line_number}: continuation event "
                        f"{continuation_event} does not match current "
                        f"event {current_event['event']}."
                    )

                add_candidate_row(
                    current_event,
                    row,
                    line_number,
                )

                parse_statistics["continuation_lines"] += 1

    if current_event is not None:
        events.append(current_event)

    return (
        events,
        reference_accession,
        parse_statistics,
    )


def evidence_tier(number_significant, support_fraction):
    """Transparent evidence classification based on independent methods."""
    if (
        number_significant >= 5
        and support_fraction >= 0.75
    ):
        return "very_high"

    if (
        number_significant >= 4
        and support_fraction >= 0.50
    ):
        return "high"

    if (
        number_significant >= 3
        and support_fraction >= 0.50
    ):
        return "moderate"

    return "low"


def finalize_event(
    event,
    sequences,
    alignment_length,
    reference_accession,
    reference_genome_length,
    args,
):
    result = {}

    result["event"] = event["event"]
    result["rdp_file_number"] = event["rdp_file_number"]

    result["reference_accession"] = reference_accession
    result["reference_genome_length"] = reference_genome_length
    result["alignment_length"] = alignment_length

    alignment_begin = parse_integer(
        event["alignment_begin_raw"]
    )
    alignment_end = parse_integer(
        event["alignment_end_raw"]
    )

    recombinant_begin = parse_integer(
        event["recombinant_begin_raw"]
    )
    recombinant_end = parse_integer(
        event["recombinant_end_raw"]
    )

    reference_begin = parse_integer(
        event["reference_begin_raw"]
    )
    reference_end = parse_integer(
        event["reference_end_raw"]
    )

    result["alignment_begin"] = alignment_begin
    result["alignment_end"] = alignment_end

    result["recombinant_begin"] = recombinant_begin
    result["recombinant_end"] = recombinant_end

    result["reference_begin"] = reference_begin
    result["reference_end"] = reference_end

    result["crosses_alignment_origin"] = (
        alignment_begin is not None
        and alignment_end is not None
        and alignment_begin > alignment_end
    )

    result["crosses_reference_origin"] = (
        reference_begin is not None
        and reference_end is not None
        and reference_begin > reference_end
    )

    (
        segment_1_start,
        segment_1_end,
        segment_2_start,
        segment_2_end,
    ) = split_circular_segment(
        reference_begin,
        reference_end,
        reference_genome_length,
    )

    result["reference_segment1_start"] = segment_1_start
    result["reference_segment1_end"] = segment_1_end
    result["reference_segment2_start"] = segment_2_start
    result["reference_segment2_end"] = segment_2_end

    fragment_length = circular_fragment_length(
        reference_begin,
        reference_end,
        reference_genome_length,
    )

    result["reference_fragment_length"] = fragment_length

    if fragment_length is not None:
        fragment_fraction = (
            fragment_length / reference_genome_length
        )
    else:
        fragment_fraction = None

    result["reference_fragment_fraction"] = fragment_fraction

    result["reference_fragment_percent"] = (
        100.0 * fragment_fraction
        if fragment_fraction is not None
        else None
    )

    recombinant_accessions = unique_preserving_order(
        event["recombinant_accessions"]
    )

    minor_parent_accessions = unique_preserving_order(
        event["minor_parent_accessions"]
    )

    major_parent_accessions = unique_preserving_order(
        event["major_parent_accessions"]
    )

    result["number_recombinants"] = len(
        recombinant_accessions
    )

    result["recombinant_accessions"] = ";".join(
        recombinant_accessions
    )

    result["recombinant_sequences_raw"] = " || ".join(
        unique_preserving_order(event["recombinant_raw"])
    )

    result["number_minor_parent_candidates"] = len(
        minor_parent_accessions
    )

    result["minor_parent_accessions"] = ";".join(
        minor_parent_accessions
    )

    result["minor_parent_candidates_raw"] = " || ".join(
        unique_preserving_order(event["minor_parent_raw"])
    )

    result["number_major_parent_candidates"] = len(
        major_parent_accessions
    )

    result["major_parent_accessions"] = ";".join(
        major_parent_accessions
    )

    result["major_parent_candidates_raw"] = " || ".join(
        unique_preserving_order(event["major_parent_raw"])
    )

    known_minor_parents = [
        accession
        for accession in minor_parent_accessions
        if accession != "Unknown"
    ]

    known_major_parents = [
        accession
        for accession in major_parent_accessions
        if accession != "Unknown"
    ]

    result["number_known_minor_parent_candidates"] = len(
        known_minor_parents
    )

    result["number_known_major_parent_candidates"] = len(
        known_major_parents
    )

    recombinant_outside_alignment = [
        accession
        for accession in recombinant_accessions
        if accession not in sequences
    ]

    parents_outside_alignment = [
        accession
        for accession in (
            known_minor_parents + known_major_parents
        )
        if accession not in sequences
    ]

    parents_outside_alignment = unique_preserving_order(
        parents_outside_alignment
    )

    result["recombinants_outside_alignment"] = ";".join(
        recombinant_outside_alignment
    )

    result["parents_outside_alignment"] = ";".join(
        parents_outside_alignment
    )

    result["all_recombinants_in_alignment"] = (
        len(recombinant_outside_alignment) == 0
    )

    result["all_known_parents_in_alignment"] = (
        len(parents_outside_alignment) == 0
    )

    independent_statuses = []
    independent_pvalues = []

    all_statuses = []
    all_pvalues = []

    for method in ALL_METHODS:
        raw_value = event["method_raw"].get(method, "")

        status, pvalue = method_status(
            raw_value,
            args.method_pvalue_cutoff,
        )

        result[f"{method}_raw"] = raw_value
        result[f"{method}_status"] = status
        result[f"{method}_pvalue"] = pvalue

        all_statuses.append(status)

        if pvalue is not None:
            all_pvalues.append(pvalue)

        if method in INDEPENDENT_METHODS:
            independent_statuses.append(status)

            if pvalue is not None:
                independent_pvalues.append(pvalue)

    number_significant = independent_statuses.count(
        "significant"
    )

    number_nonsignificant = independent_statuses.count(
        "nonsignificant"
    )

    number_tested = (
        number_significant + number_nonsignificant
    )

    number_not_reported = len(
        INDEPENDENT_METHODS
    ) - number_tested

    support_fraction = (
        number_significant / number_tested
        if number_tested > 0
        else 0.0
    )

    result["number_significant_methods"] = (
        number_significant
    )

    result["number_nonsignificant_methods"] = (
        number_nonsignificant
    )

    result["number_tested_methods"] = number_tested
    result["number_not_reported_methods"] = (
        number_not_reported
    )

    result["support_fraction"] = support_fraction

    result["number_significant_methods_all9"] = (
        all_statuses.count("significant")
    )

    result["number_tested_methods_all9"] = (
        all_statuses.count("significant")
        + all_statuses.count("nonsignificant")
    )

    result["minimum_pvalue"] = (
        min(independent_pvalues)
        if independent_pvalues
        else None
    )

    positive_pvalues = [
        value
        for value in independent_pvalues
        if value > 0
    ]

    if any(value == 0 for value in independent_pvalues):
        result["maximum_minus_log10_p"] = "inf"
    elif positive_pvalues:
        result["maximum_minus_log10_p"] = max(
            -math.log10(value)
            for value in positive_pvalues
        )
    else:
        result["maximum_minus_log10_p"] = None

    result["evidence_tier"] = evidence_tier(
        number_significant,
        support_fraction,
    )

    result["contains_alternative_process_marker"] = (
        event["contains_alternative_process_marker"]
    )

    result["contains_undetermined_breakpoint_marker"] = (
        event["contains_undetermined_breakpoint_marker"]
    )

    result["contains_recombinant_caret_marker"] = (
        event["contains_recombinant_caret_marker"]
    )

    result["contains_parental_marker"] = (
        event["contains_parental_marker"]
    )

    result["contains_unknown_minor_parent"] = (
        event["contains_unknown_minor_parent"]
    )

    result["contains_unknown_major_parent"] = (
        event["contains_unknown_major_parent"]
    )

    has_complete_reference_coordinates = (
        reference_begin is not None
        and reference_end is not None
        and reference_genome_length is not None
    )

    result["has_complete_reference_coordinates"] = (
        has_complete_reference_coordinates
    )

    mapping_exclusion_reasons = []

    if not has_complete_reference_coordinates:
        mapping_exclusion_reasons.append(
            "incomplete_reference_coordinates"
        )

    if event["contains_undetermined_breakpoint_marker"]:
        mapping_exclusion_reasons.append(
            "undetermined_breakpoint"
        )

    if event["contains_recombinant_caret_marker"]:
        mapping_exclusion_reasons.append(
            "uncertain_recombinant_identity"
        )

    if event["contains_alternative_process_marker"]:
        mapping_exclusion_reasons.append(
            "possible_non_recombination_process"
        )

    result["eligible_for_breakpoint_comparison"] = (
        len(mapping_exclusion_reasons) == 0
    )

    result["mapping_exclusion_reason"] = ";".join(
        mapping_exclusion_reasons
    )

    result["criterion_min_methods"] = args.min_methods
    result["criterion_min_tested_methods"] = (
        args.min_tested_methods
    )
    result["criterion_support_fraction"] = (
        args.min_support_fraction
    )
    result["criterion_max_missing_methods"] = (
        args.max_missing_methods
    )
    result["criterion_method_pvalue_cutoff"] = (
        args.method_pvalue_cutoff
    )
    result["criterion_min_fragment_length"] = (
        args.min_fragment_length
    )
    result["criterion_internal_3seq_excluded"] = True
    result["criterion_allow_caret"] = args.allow_caret
    result["criterion_allow_tilde"] = args.allow_tilde
    result["criterion_require_known_parents"] = (
        args.require_known_parents
    )

    failed_criteria = []

    if number_significant < args.min_methods:
        failed_criteria.append("min_methods")

    if number_tested < args.min_tested_methods:
        failed_criteria.append("min_tested_methods")

    if support_fraction < args.min_support_fraction:
        failed_criteria.append("support_fraction")

    if number_not_reported > args.max_missing_methods:
        failed_criteria.append("max_missing_methods")

    if (
        fragment_length is None
        or fragment_length < args.min_fragment_length
    ):
        failed_criteria.append("min_fragment_length")

    if (
        event["contains_recombinant_caret_marker"]
        and not args.allow_caret
    ):
        failed_criteria.append("caret_marker")

    if (
        event["contains_alternative_process_marker"]
        and not args.allow_tilde
    ):
        failed_criteria.append("tilde_marker")

    if args.require_known_parents:
        if not known_minor_parents:
            failed_criteria.append("known_minor_parent")

        if not known_major_parents:
            failed_criteria.append("known_major_parent")

    if recombinant_outside_alignment:
        failed_criteria.append(
            "recombinant_outside_alignment"
        )

    result["robust_event"] = len(failed_criteria) == 0

    result["failed_robustness_criteria"] = ";".join(
        failed_criteria
    )

    return result


def event_fieldnames():
    fields = [
        "event",
        "rdp_file_number",
        "reference_accession",
        "reference_genome_length",
        "alignment_length",
        "alignment_begin",
        "alignment_end",
        "crosses_alignment_origin",
        "recombinant_begin",
        "recombinant_end",
        "reference_begin",
        "reference_end",
        "crosses_reference_origin",
        "reference_segment1_start",
        "reference_segment1_end",
        "reference_segment2_start",
        "reference_segment2_end",
        "reference_fragment_length",
        "reference_fragment_fraction",
        "reference_fragment_percent",
        "number_recombinants",
        "recombinant_accessions",
        "recombinant_sequences_raw",
        "number_minor_parent_candidates",
        "number_known_minor_parent_candidates",
        "minor_parent_accessions",
        "minor_parent_candidates_raw",
        "number_major_parent_candidates",
        "number_known_major_parent_candidates",
        "major_parent_accessions",
        "major_parent_candidates_raw",
        "recombinants_outside_alignment",
        "parents_outside_alignment",
        "all_recombinants_in_alignment",
        "all_known_parents_in_alignment",
    ]

    for method in ALL_METHODS:
        fields.extend(
            [
                f"{method}_raw",
                f"{method}_status",
                f"{method}_pvalue",
            ]
        )

    fields.extend(
        [
            "number_significant_methods",
            "number_nonsignificant_methods",
            "number_tested_methods",
            "number_not_reported_methods",
            "support_fraction",
            "number_significant_methods_all9",
            "number_tested_methods_all9",
            "minimum_pvalue",
            "maximum_minus_log10_p",
            "evidence_tier",
            "contains_alternative_process_marker",
            "contains_undetermined_breakpoint_marker",
            "contains_recombinant_caret_marker",
            "contains_parental_marker",
            "contains_unknown_minor_parent",
            "contains_unknown_major_parent",
            "has_complete_reference_coordinates",
            "eligible_for_breakpoint_comparison",
            "mapping_exclusion_reason",
            "criterion_min_methods",
            "criterion_min_tested_methods",
            "criterion_support_fraction",
            "criterion_max_missing_methods",
            "criterion_method_pvalue_cutoff",
            "criterion_min_fragment_length",
            "criterion_internal_3seq_excluded",
            "criterion_allow_caret",
            "criterion_allow_tilde",
            "criterion_require_known_parents",
            "robust_event",
            "failed_robustness_criteria",
        ]
    )

    return fields


def write_event_table(filename, events, only_robust=False):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = event_fieldnames()

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()

        for event in events:
            if only_robust and not event["robust_event"]:
                continue

            writer.writerow(event)


def write_candidate_table(filename, raw_events):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "event",
        "rdp_file_number",
        "source_line",
        "recombinant_raw",
        "recombinant_accession",
        "minor_parent_raw",
        "minor_parent_accession",
        "major_parent_raw",
        "major_parent_accession",
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for event in raw_events:
            for candidate in event["candidate_rows"]:
                row = {
                    "event": event["event"],
                    "rdp_file_number": event["rdp_file_number"],
                    **candidate,
                }

                writer.writerow(row)


def write_summary(
    filename,
    events,
    reference_accession,
    reference_genome_length,
    alignment_length,
    parse_statistics,
    args,
):
    robust_events = [
        event
        for event in events
        if event["robust_event"]
    ]

    summary = {
        "input_report": str(Path(args.input).resolve()),
        "alignment": str(Path(args.alignment).resolve()),
        "reference_accession": reference_accession,
        "reference_genome_length": reference_genome_length,
        "alignment_length": alignment_length,
        "number_alignment_sequences": None,
        "number_events_total": len(events),
        "number_events_robust": len(robust_events),
        "number_events_not_robust": (
            len(events) - len(robust_events)
        ),
        "number_events_eligible_for_breakpoint_comparison": sum(
            event["eligible_for_breakpoint_comparison"]
            for event in events
        ),
        "number_events_with_caret_marker": sum(
            event["contains_recombinant_caret_marker"]
            for event in events
        ),
        "number_events_with_tilde_marker": sum(
            event["contains_alternative_process_marker"]
            for event in events
        ),
        "number_events_with_undetermined_breakpoints": sum(
            event["contains_undetermined_breakpoint_marker"]
            for event in events
        ),
        "number_events_with_unknown_minor_parent": sum(
            event["contains_unknown_minor_parent"]
            for event in events
        ),
        "number_events_with_unknown_major_parent": sum(
            event["contains_unknown_major_parent"]
            for event in events
        ),
        "independent_methods": INDEPENDENT_METHODS,
        "internal_3seq_excluded": True,
        "criteria": {
            "min_methods": args.min_methods,
            "min_tested_methods": args.min_tested_methods,
            "min_support_fraction": args.min_support_fraction,
            "max_missing_methods": args.max_missing_methods,
            "method_pvalue_cutoff": args.method_pvalue_cutoff,
            "min_fragment_length": args.min_fragment_length,
            "allow_caret": args.allow_caret,
            "allow_tilde": args.allow_tilde,
            "require_known_parents": args.require_known_parents,
        },
        "parse_statistics": parse_statistics,
    }

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            summary,
            handle,
            indent=2,
            ensure_ascii=False,
        )


def validate_arguments(args):
    if args.min_methods < 1:
        raise ValueError("--min-methods must be positive.")

    if not 1 <= args.min_tested_methods <= len(
        INDEPENDENT_METHODS
    ):
        raise ValueError(
            "--min-tested-methods must be between 1 and 8."
        )

    if not 0 <= args.min_support_fraction <= 1:
        raise ValueError(
            "--min-support-fraction must be between 0 and 1."
        )

    if not 0 <= args.max_missing_methods <= len(
        INDEPENDENT_METHODS
    ):
        raise ValueError(
            "--max-missing-methods must be between 0 and 8."
        )

    if not 0 < args.method_pvalue_cutoff <= 1:
        raise ValueError(
            "--method-pvalue-cutoff must be in ]0, 1]."
        )

    if args.min_fragment_length < 1:
        raise ValueError(
            "--min-fragment-length must be positive."
        )


def main():
    args = parse_arguments()
    validate_arguments(args)

    sequences, alignment_length = read_alignment(
        args.alignment
    )

    (
        raw_events,
        reference_accession,
        parse_statistics,
    ) = parse_rdp_report(
        args.input,
        reference_override=args.reference_accession,
    )

    reference_genome_length = (
        reference_length_from_alignment(
            sequences,
            reference_accession,
        )
    )

    events = [
        finalize_event(
            event=event,
            sequences=sequences,
            alignment_length=alignment_length,
            reference_accession=reference_accession,
            reference_genome_length=reference_genome_length,
            args=args,
        )
        for event in raw_events
    ]

    write_event_table(
        args.events_all,
        events,
        only_robust=False,
    )

    write_event_table(
        args.events_robust,
        events,
        only_robust=True,
    )

    write_candidate_table(
        args.candidates,
        raw_events,
    )

    write_summary(
        args.summary,
        events,
        reference_accession,
        reference_genome_length,
        alignment_length,
        parse_statistics,
        args,
    )

    number_robust = sum(
        event["robust_event"]
        for event in events
    )

    print(f"Reference accession       : {reference_accession}")
    print(f"Reference genome length   : {reference_genome_length}")
    print(f"Alignment length          : {alignment_length}")
    print(f"Alignment sequences       : {len(sequences)}")
    print(f"Main RDP events parsed    : {len(events)}")
    print(f"Robust independent events : {number_robust}")
    print(
        "Breakpoint-comparable    : "
        f"{sum(event['eligible_for_breakpoint_comparison'] for event in events)}"
    )
    print(f"Continuation lines        : {parse_statistics['continuation_lines']}")
    print(f"Blank lines ignored       : {parse_statistics['blank_lines']}")


if __name__ == "__main__":
    main()
