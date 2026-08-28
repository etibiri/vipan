#!/usr/bin/env python3

import argparse
import csv
import math
import re
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path


INTEGER_PATTERN = re.compile(r"-?\d+")
THRESHOLD_PATTERN = re.compile(
    r"Need a p-value of\s+([0-9.eE+-]+)",
    flags=re.IGNORECASE,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Collapse redundant 3SEQ recombinant triplets into "
            "non-redundant breakpoint events."
        )
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--run-log", required=True)
    parser.add_argument("--events", required=True)
    parser.add_argument("--summary", required=True)

    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Minimum recombinant segment length.",
    )

    parser.add_argument(
        "--breakpoint-bin",
        type=int,
        default=25,
        help="Breakpoint clustering bin size in nucleotides.",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=250000,
        help="Number of input lines processed before each database flush.",
    )

    parser.add_argument(
        "--tmpdir",
        default="/tmp",
        help="Temporary directory used for the SQLite database.",
    )

    return parser.parse_args()


def get_alignment_length(fasta_file):
    """Return the aligned length of the first FASTA sequence."""
    sequence = []
    sequence_started = False

    with open(fasta_file, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence_started:
                    break

                sequence_started = True
                continue

            if sequence_started:
                sequence.append(line)

    if not sequence:
        raise ValueError(
            f"No sequence found in alignment: {fasta_file}"
        )

    return len("".join(sequence))


def get_significance_threshold(log_file):
    """Extract the multiple-testing threshold from the 3SEQ log."""
    log_path = Path(log_file)

    if not log_path.exists():
        return None

    text = log_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    match = THRESHOLD_PATTERN.search(text)

    if match is None:
        return None

    return float(match.group(1))


def parse_integer(value, default=0):
    match = INTEGER_PATTERN.search(str(value))

    if match is None:
        return default

    return int(match.group())


def parse_float(value, default=math.inf):
    text = str(value).strip().lower()

    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def normalize_position(position, genome_length):
    """Normalize a coordinate to the interval 1..genome_length."""
    return ((int(position) - 1) % genome_length) + 1


def circular_midpoint(start, end, genome_length):
    """
    Calculate the midpoint of a potentially origin-crossing interval.
    """
    start = normalize_position(start, genome_length)
    end = normalize_position(end, genome_length)

    forward_distance = (end - start) % genome_length
    midpoint = start + round(forward_distance / 2)

    return normalize_position(midpoint, genome_length)


def parse_breakpoints(text, genome_length):
    """
    Parse 3SEQ breakpoint coordinates.

    Two numbers are interpreted as two exact breakpoint positions.
    Four or more numbers are interpreted as two breakpoint ranges:
    start1-end1 and start2-end2.
    """
    values = [
        int(value)
        for value in INTEGER_PATTERN.findall(text)
    ]

    if len(values) == 2:
        breakpoint_1 = normalize_position(
            values[0],
            genome_length,
        )
        breakpoint_2 = normalize_position(
            values[1],
            genome_length,
        )

    elif len(values) >= 4:
        breakpoint_1 = circular_midpoint(
            values[0],
            values[1],
            genome_length,
        )
        breakpoint_2 = circular_midpoint(
            values[2],
            values[3],
            genome_length,
        )

    else:
        return None

    return breakpoint_1, breakpoint_2


def breakpoint_bin(position, genome_length, requested_bin_size):
    """
    Assign a position to a circular breakpoint bin.

    The first and last genomic bins are adjacent, which is important
    for circular viral genomes.
    """
    number_of_bins = max(
        1,
        round(genome_length / requested_bin_size),
    )

    scaled_position = (
        (position - 1)
        * number_of_bins
        / genome_length
    )

    bin_identifier = int(
        math.floor(scaled_position + 0.5)
    ) % number_of_bins

    return bin_identifier, number_of_bins


def bin_center(bin_identifier, genome_length, number_of_bins):
    center = round(
        bin_identifier
        * genome_length
        / number_of_bins
    ) + 1

    return normalize_position(center, genome_length)


def is_hs_approximation(value):
    value = str(value).strip().upper()

    return value in {
        "Y",
        "YES",
        "TRUE",
        "1",
        "HS",
        "*",
    }


def create_database(connection):
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        PRAGMA cache_size = -200000;

        CREATE TABLE event_counts (
            child TEXT NOT NULL,
            bin1 INTEGER NOT NULL,
            bin2 INTEGER NOT NULL,
            supporting_triplets INTEGER NOT NULL,
            PRIMARY KEY (child, bin1, bin2)
        );

        CREATE TABLE representatives (
            child TEXT NOT NULL,
            bin1 INTEGER NOT NULL,
            bin2 INTEGER NOT NULL,

            parent1 TEXT NOT NULL,
            parent2 TEXT NOT NULL,

            breakpoint1 INTEGER NOT NULL,
            breakpoint2 INTEGER NOT NULL,
            breakpoint_raw TEXT NOT NULL,

            m INTEGER NOT NULL,
            n INTEGER NOT NULL,
            k INTEGER NOT NULL,

            p_raw TEXT NOT NULL,
            p_numeric REAL NOT NULL,

            hs_raw TEXT NOT NULL,
            hs_rank INTEGER NOT NULL,

            log_p_raw TEXT NOT NULL,
            log_p_numeric REAL NOT NULL,

            ds_p1 TEXT NOT NULL,
            ds_p2 TEXT NOT NULL,

            min_recombinant_length INTEGER NOT NULL,

            PRIMARY KEY (child, bin1, bin2)
        );
        """
    )


REPRESENTATIVE_UPSERT = """
INSERT INTO representatives (
    child,
    bin1,
    bin2,
    parent1,
    parent2,
    breakpoint1,
    breakpoint2,
    breakpoint_raw,
    m,
    n,
    k,
    p_raw,
    p_numeric,
    hs_raw,
    hs_rank,
    log_p_raw,
    log_p_numeric,
    ds_p1,
    ds_p2,
    min_recombinant_length
)
VALUES (
    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
)
ON CONFLICT (child, bin1, bin2)
DO UPDATE SET
    parent1 = excluded.parent1,
    parent2 = excluded.parent2,
    breakpoint1 = excluded.breakpoint1,
    breakpoint2 = excluded.breakpoint2,
    breakpoint_raw = excluded.breakpoint_raw,
    m = excluded.m,
    n = excluded.n,
    k = excluded.k,
    p_raw = excluded.p_raw,
    p_numeric = excluded.p_numeric,
    hs_raw = excluded.hs_raw,
    hs_rank = excluded.hs_rank,
    log_p_raw = excluded.log_p_raw,
    log_p_numeric = excluded.log_p_numeric,
    ds_p1 = excluded.ds_p1,
    ds_p2 = excluded.ds_p2,
    min_recombinant_length =
        excluded.min_recombinant_length
WHERE
       excluded.p_numeric < representatives.p_numeric

    OR (
        excluded.p_numeric = representatives.p_numeric
        AND excluded.hs_rank < representatives.hs_rank
    )

    OR (
        excluded.p_numeric = representatives.p_numeric
        AND excluded.hs_rank = representatives.hs_rank
        AND excluded.log_p_numeric
            < representatives.log_p_numeric
    )

    OR (
        excluded.p_numeric = representatives.p_numeric
        AND excluded.hs_rank = representatives.hs_rank
        AND excluded.log_p_numeric
            = representatives.log_p_numeric
        AND excluded.k > representatives.k
    )

    OR (
        excluded.p_numeric = representatives.p_numeric
        AND excluded.hs_rank = representatives.hs_rank
        AND excluded.log_p_numeric
            = representatives.log_p_numeric
        AND excluded.k = representatives.k
        AND excluded.min_recombinant_length
            > representatives.min_recombinant_length
    );
"""


COUNT_UPSERT = """
INSERT INTO event_counts (
    child,
    bin1,
    bin2,
    supporting_triplets
)
VALUES (?, ?, ?, ?)
ON CONFLICT (child, bin1, bin2)
DO UPDATE SET
    supporting_triplets =
        event_counts.supporting_triplets
        + excluded.supporting_triplets;
"""


def flush_buffers(
    connection,
    count_buffer,
    representative_buffer,
):
    if not count_buffer:
        return

    count_rows = [
        (
            child,
            bin1,
            bin2,
            count,
        )
        for (child, bin1, bin2), count
        in count_buffer.items()
    ]

    representative_rows = [
        record
        for _, record in representative_buffer.values()
    ]

    connection.executemany(
        COUNT_UPSERT,
        count_rows,
    )

    connection.executemany(
        REPRESENTATIVE_UPSERT,
        representative_rows,
    )

    connection.commit()

    count_buffer.clear()
    representative_buffer.clear()


def process_input(
    input_file,
    connection,
    genome_length,
    requested_bin_size,
    min_length,
    significance_threshold,
    chunk_size,
):
    count_buffer = defaultdict(int)
    representative_buffer = {}

    statistics = {
        "total_lines": 0,
        "retained_lines": 0,
        "short_segments": 0,
        "non_significant": 0,
        "malformed_lines": 0,
        "unparsed_breakpoints": 0,
    }

    number_of_bins = max(
        1,
        round(genome_length / requested_bin_size),
    )

    with open(
        input_file,
        encoding="utf-8",
        errors="replace",
    ) as handle:
        header = handle.readline().rstrip("\n").split("\t")

        if len(header) < 13:
            raise ValueError(
                "Unexpected 3SEQ header: fewer than 13 columns."
            )

        for line in handle:
            statistics["total_lines"] += 1

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) < 13:
                statistics["malformed_lines"] += 1
                continue

            parent1 = fields[0].strip()
            parent2 = fields[1].strip()
            child = fields[2].strip()

            minimum_length = parse_integer(
                fields[11],
                default=-1,
            )

            if minimum_length < min_length:
                statistics["short_segments"] += 1
                continue

            p_raw = fields[6].strip()
            p_numeric = parse_float(p_raw)

            if not math.isfinite(p_numeric):
                statistics["malformed_lines"] += 1
                continue

            if (
                significance_threshold is not None
                and p_numeric > significance_threshold
            ):
                statistics["non_significant"] += 1
                continue

            breakpoint_raw = " ".join(
                fields[12:]
            ).strip()

            parsed_breakpoints = parse_breakpoints(
                breakpoint_raw,
                genome_length,
            )

            if parsed_breakpoints is None:
                statistics["unparsed_breakpoints"] += 1
                continue

            breakpoint_1, breakpoint_2 = parsed_breakpoints

            bin_1, _ = breakpoint_bin(
                breakpoint_1,
                genome_length,
                requested_bin_size,
            )

            bin_2, _ = breakpoint_bin(
                breakpoint_2,
                genome_length,
                requested_bin_size,
            )

            # Treat the pair of breakpoints as unordered.
            if (bin_1, breakpoint_1) > (bin_2, breakpoint_2):
                bin_1, bin_2 = bin_2, bin_1
                breakpoint_1, breakpoint_2 = (
                    breakpoint_2,
                    breakpoint_1,
                )

            event_key = (
                child,
                bin_1,
                bin_2,
            )

            m_value = parse_integer(fields[3])
            n_value = parse_integer(fields[4])
            k_value = parse_integer(fields[5])

            hs_raw = fields[7].strip()
            hs_rank = int(
                is_hs_approximation(hs_raw)
            )

            log_p_raw = fields[8].strip()
            log_p_numeric = parse_float(
                log_p_raw,
                default=math.inf,
            )

            ds_p1 = fields[9].strip()
            ds_p2 = fields[10].strip()

            ranking = (
                p_numeric,
                hs_rank,
                log_p_numeric,
                -k_value,
                -minimum_length,
            )

            record = (
                child,
                bin_1,
                bin_2,
                parent1,
                parent2,
                breakpoint_1,
                breakpoint_2,
                breakpoint_raw,
                m_value,
                n_value,
                k_value,
                p_raw,
                p_numeric,
                hs_raw,
                hs_rank,
                log_p_raw,
                log_p_numeric,
                ds_p1,
                ds_p2,
                minimum_length,
            )

            count_buffer[event_key] += 1

            previous = representative_buffer.get(
                event_key
            )

            if previous is None or ranking < previous[0]:
                representative_buffer[event_key] = (
                    ranking,
                    record,
                )

            statistics["retained_lines"] += 1

            if (
                statistics["total_lines"]
                % chunk_size
                == 0
            ):
                flush_buffers(
                    connection,
                    count_buffer,
                    representative_buffer,
                )

    flush_buffers(
        connection,
        count_buffer,
        representative_buffer,
    )

    return statistics, number_of_bins


def write_results(
    connection,
    events_file,
    summary_file,
    genome_length,
    requested_bin_size,
    number_of_bins,
    significance_threshold,
):
    events_path = Path(events_file)
    summary_path = Path(summary_file)

    events_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summaries = {}

    query = """
    SELECT
        representatives.*,
        event_counts.supporting_triplets
    FROM representatives
    INNER JOIN event_counts
        USING (child, bin1, bin2)
    ORDER BY
        child,
        bin1,
        bin2;
    """

    connection.row_factory = sqlite3.Row

    threshold_text = (
        f"{significance_threshold:.8e}"
        if significance_threshold is not None
        else "NA"
    )

    event_count = 0

    with events_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "Event_ID",
                "Recombinant",
                "Best_parent_1",
                "Best_parent_2",
                "Breakpoint_1",
                "Breakpoint_2",
                "Breakpoint_range_raw",
                "Breakpoint_cluster_1",
                "Breakpoint_cluster_2",
                "Breakpoint_bin_size",
                "Min_recombinant_length",
                "m",
                "n",
                "k",
                "P_value",
                "HS_approximation",
                "log_P",
                "DS_P_1",
                "DS_P_2",
                "Supporting_triplets",
                "Multiple_testing_threshold",
            ]
        )

        for row in connection.execute(query):
            event_count += 1
            event_id = f"3SEQ_E{event_count:06d}"

            center_1 = bin_center(
                row["bin1"],
                genome_length,
                number_of_bins,
            )

            center_2 = bin_center(
                row["bin2"],
                genome_length,
                number_of_bins,
            )

            writer.writerow(
                [
                    event_id,
                    row["child"],
                    row["parent1"],
                    row["parent2"],
                    row["breakpoint1"],
                    row["breakpoint2"],
                    row["breakpoint_raw"],
                    center_1,
                    center_2,
                    requested_bin_size,
                    row["min_recombinant_length"],
                    row["m"],
                    row["n"],
                    row["k"],
                    row["p_raw"],
                    row["hs_raw"],
                    row["log_p_raw"],
                    row["ds_p1"],
                    row["ds_p2"],
                    row["supporting_triplets"],
                    threshold_text,
                ]
            )

            summary = summaries.setdefault(
                row["child"],
                {
                    "events": 0,
                    "triplets": 0,
                    "best_rank": None,
                    "best_p": None,
                    "best_log_p": None,
                    "maximum_length": 0,
                },
            )

            summary["events"] += 1
            summary["triplets"] += row[
                "supporting_triplets"
            ]

            rank = (
                row["p_numeric"],
                row["hs_rank"],
                row["log_p_numeric"],
                -row["k"],
                -row["min_recombinant_length"],
            )

            if (
                summary["best_rank"] is None
                or rank < summary["best_rank"]
            ):
                summary["best_rank"] = rank
                summary["best_p"] = row["p_raw"]
                summary["best_log_p"] = row["log_p_raw"]

            summary["maximum_length"] = max(
                summary["maximum_length"],
                row["min_recombinant_length"],
            )

    with summary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "Recombinant",
                "Nonredundant_events",
                "Supporting_triplets",
                "Best_P_value",
                "Best_log_P",
                "Maximum_min_recombinant_length",
            ]
        )

        for recombinant, summary in summaries.items():
            writer.writerow(
                [
                    recombinant,
                    summary["events"],
                    summary["triplets"],
                    summary["best_p"],
                    summary["best_log_p"],
                    summary["maximum_length"],
                ]
            )

    return event_count, len(summaries)


def main():
    args = parse_arguments()

    if args.min_length < 1:
        raise ValueError("--min-length must be positive.")

    if args.breakpoint_bin < 1:
        raise ValueError("--breakpoint-bin must be positive.")

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be positive.")

    genome_length = get_alignment_length(
        args.alignment
    )

    threshold = get_significance_threshold(
        args.run_log
    )

    print(f"Alignment length: {genome_length}")
    print(
        "Multiple-testing threshold: "
        f"{threshold if threshold is not None else 'not found'}"
    )
    print(
        f"Minimum recombinant length: {args.min_length}"
    )
    print(
        f"Breakpoint bin size: {args.breakpoint_bin}"
    )

    with tempfile.TemporaryDirectory(
        prefix="collapse_3seq_",
        dir=args.tmpdir,
    ) as temporary_directory:
        database_file = (
            Path(temporary_directory)
            / "events.sqlite"
        )

        connection = sqlite3.connect(database_file)

        try:
            create_database(connection)

            statistics, number_of_bins = process_input(
                input_file=args.input,
                connection=connection,
                genome_length=genome_length,
                requested_bin_size=args.breakpoint_bin,
                min_length=args.min_length,
                significance_threshold=threshold,
                chunk_size=args.chunk_size,
            )

            event_count, recombinant_count = write_results(
                connection=connection,
                events_file=args.events,
                summary_file=args.summary,
                genome_length=genome_length,
                requested_bin_size=args.breakpoint_bin,
                number_of_bins=number_of_bins,
                significance_threshold=threshold,
            )

        finally:
            connection.close()

    print(f"Input triplets: {statistics['total_lines']}")
    print(f"Retained triplets: {statistics['retained_lines']}")
    print(f"Short segments removed: {statistics['short_segments']}")
    print(
        "Non-significant triplets removed: "
        f"{statistics['non_significant']}"
    )
    print(
        "Malformed lines removed: "
        f"{statistics['malformed_lines']}"
    )
    print(
        "Unparsed breakpoint lines removed: "
        f"{statistics['unparsed_breakpoints']}"
    )
    print(f"Non-redundant events: {event_count}")
    print(f"Recombinant sequences: {recombinant_count}")


if __name__ == "__main__":
    main()
