"""Small dependency-free FASTA helpers shared by workflow utility scripts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List


@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def read_fasta(path: Path) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header = None
    sequence_parts: List[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    identifier = header.split()[0]
                    records.append(FastaRecord(identifier, header, "".join(sequence_parts)))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header at line {line_number} in {path}")
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"sequence before first FASTA header at line {line_number} in {path}")
                sequence_parts.append("".join(line.split()))
    if header is not None:
        records.append(FastaRecord(header.split()[0], header, "".join(sequence_parts)))
    return records


def write_fasta(records: Iterable[FastaRecord], path: Path, width: int = 80) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            count += 1
            handle.write(f">{record.id}\n")
            sequence = record.sequence
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")
    return count


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return sequence.translate(translation)[::-1]
