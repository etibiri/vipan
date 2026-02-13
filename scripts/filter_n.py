#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import os

def filter_sequences(input_fasta, output_fasta, max_n_pct):
    kept = 0
    total = 0
    
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq = str(record.seq).upper()
            n_count = seq.count('N')
            n_pct = (n_count / len(seq)) * 100
            
            if n_pct <= max_n_pct:
                SeqIO.write(record, out_handle, "fasta")
                kept += 1
    
    return total, kept

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences with too many Ns")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max_n_pct", type=float, default=1.0)
    args = parser.parse_args()

    total, kept = filter_sequences(args.input, args.output, args.max_n_pct)
    print(f"Total sequences: {total}")
    print(f"Kept (N <= {args.max_n_pct}%): {kept}")
    print(f"Removed: {total - kept}")
