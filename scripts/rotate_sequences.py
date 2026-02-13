#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def rotate_to_ori(input_fasta, output_fasta, motif="TAATATTAC"):
    rotated_count = 0
    total = 0
    
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq_str = str(record.seq).upper()
            
            # Find the motif in the sequence
            pos = seq_str.find(motif)
            
            if pos != -1:
                # Rotate: everything from pos to end + everything from start to pos
                new_seq = seq_str[pos:] + seq_str[:pos]
                new_record = SeqRecord(
                    record.seq.__class__(new_seq),
                    id=record.id,
                    description=record.description
                )
                SeqIO.write(new_record, out_handle, "fasta")
                rotated_count += 1
            else:
                # If motif not found, keep as is (or you could log a warning)
                SeqIO.write(record, out_handle, "fasta")
    
    return total, rotated_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rotate circular genomes to a specific motif")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--motif", default="TAATATTAC")
    args = parser.parse_args()

    total, rotated = rotate_to_ori(args.input, args.output, args.motif)
    print(f"Total sequences: {total}")
    print(f"Rotated to {args.motif}: {rotated}")
