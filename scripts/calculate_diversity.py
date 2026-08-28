#!/usr/bin/env python3

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO


VALID_BASES = frozenset("ACGT")


def calculate_consensus(alignment, threshold=0.5):
    """
    Generate a nucleotide consensus sequence from a multiple alignment.

    Gaps and ambiguous characters are ignored. A position is represented
    by N when:
      - no valid A/C/G/T nucleotide is available;
      - the most frequent nucleotide does not reach the threshold;
      - several nucleotides are tied for the highest frequency.
    """
    if not 0 < threshold <= 1:
        raise ValueError("Consensus threshold must be between 0 and 1.")

    consensus = []
    alignment_length = alignment.get_alignment_length()

    for position in range(alignment_length):
        column = alignment[:, position].upper()
        bases = [base for base in column if base in VALID_BASES]

        if not bases:
            consensus.append("N")
            continue

        counts = Counter(bases)
        maximum_count = max(counts.values())

        most_frequent_bases = [
            base
            for base, count in counts.items()
            if count == maximum_count
        ]

        frequency = maximum_count / len(bases)

        if len(most_frequent_bases) == 1 and frequency >= threshold:
            consensus.append(most_frequent_bases[0])
        else:
            consensus.append("N")

    return "".join(consensus)


def calculate_site_pi(column):
    """
    Calculate unbiased nucleotide diversity for one alignment column.

    Gaps and ambiguous nucleotides are treated as missing data.
    NaN is returned when fewer than two valid nucleotides are available.
    """
    bases = [base for base in column.upper() if base in VALID_BASES]
    number_of_bases = len(bases)

    if number_of_bases < 2:
        return np.nan

    counts = np.asarray(
        list(Counter(bases).values()),
        dtype=float,
    )

    frequencies = counts / number_of_bases

    # Unbiased nucleotide diversity:
    # pi = n/(n-1) * (1 - sum(p_i^2))
    pi = (
        number_of_bases
        / (number_of_bases - 1)
        * (1.0 - np.sum(frequencies**2))
    )

    return float(pi)


def calculate_pi(alignment, window_size=50):
    """
    Calculate nucleotide diversity along an alignment using sliding windows.

    The value for each window is the mean of the callable sites. Sites with
    fewer than two valid A/C/G/T nucleotides are excluded from the mean.
    """
    number_of_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    if number_of_sequences < 2:
        raise ValueError(
            "At least two aligned sequences are required to calculate pi."
        )

    if window_size < 1:
        raise ValueError("Window size must be greater than or equal to 1.")

    if window_size > alignment_length:
        raise ValueError(
            f"Window size ({window_size}) exceeds alignment length "
            f"({alignment_length})."
        )

    # Calculate pi once for every alignment position.
    site_pi_values = np.asarray(
        [
            calculate_site_pi(alignment[:, position])
            for position in range(alignment_length)
        ],
        dtype=float,
    )

    window_pi_values = []
    callable_site_counts = []

    # +1 ensures inclusion of the final complete window.
    number_of_windows = alignment_length - window_size + 1

    for start in range(number_of_windows):
        end = start + window_size
        window_values = site_pi_values[start:end]

        callable_mask = np.isfinite(window_values)
        callable_sites = int(np.sum(callable_mask))

        if callable_sites == 0:
            mean_pi = np.nan
        else:
            mean_pi = float(np.mean(window_values[callable_mask]))

        window_pi_values.append(mean_pi)
        callable_site_counts.append(callable_sites)

    return window_pi_values, callable_site_counts


def write_consensus(consensus, output_file, line_width=80):
    """Write the consensus sequence in FASTA format."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(">Pangenome_Consensus\n")

        for start in range(0, len(consensus), line_width):
            handle.write(consensus[start:start + line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a consensus sequence and calculate sliding-window "
            "nucleotide diversity from a multiple FASTA alignment."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input multiple sequence alignment in FASTA format.",
    )

    parser.add_argument(
        "--consensus",
        required=True,
        help="Output consensus sequence in FASTA format.",
    )

    parser.add_argument(
        "--plot_data",
        required=True,
        help="Output nucleotide diversity table in TSV format.",
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=50,
        help="Sliding-window size in alignment positions (default: 50).",
    )

    parser.add_argument(
        "--consensus-threshold",
        type=float,
        default=0.5,
        help=(
            "Minimum frequency required to call a consensus nucleotide "
            "(default: 0.5)."
        ),
    )

    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")

    if len(alignment) == 0:
        raise ValueError("The input alignment contains no sequences.")

    alignment_length = alignment.get_alignment_length()

    if alignment_length == 0:
        raise ValueError("The input alignment has a length of zero.")

    print(f"Number of sequences: {len(alignment)}")
    print(f"Alignment length: {alignment_length}")
    print(f"Sliding-window size: {args.window_size}")

    # 1. Generate the consensus sequence.
    consensus = calculate_consensus(
        alignment,
        threshold=args.consensus_threshold,
    )

    write_consensus(
        consensus,
        args.consensus,
    )

    # 2. Calculate sliding-window nucleotide diversity.
    pi_values, callable_sites = calculate_pi(
        alignment,
        window_size=args.window_size,
    )

    # Position corresponds to the 1-based beginning of each window.
    diversity_table = pd.DataFrame(
        {
            "Position": np.arange(1, len(pi_values) + 1),
            "Window_end": np.arange(
                args.window_size,
                args.window_size + len(pi_values),
            ),
            "Pi": pi_values,
            "Callable_sites": callable_sites,
        }
    )

    output_path = Path(args.plot_data)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    diversity_table.to_csv(
        output_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(f"Consensus written to: {args.consensus}")
    print(f"Nucleotide diversity written to: {args.plot_data}")


if __name__ == "__main__":
    main()
