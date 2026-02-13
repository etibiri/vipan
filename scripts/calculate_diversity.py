#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd
import numpy as np

def calculate_pi(alignment, window_size=50):
    """Calculate Nucleotide Diversity (pi) along the alignment."""
    num_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()
    pi_values = []

    for i in range(0, aln_len - window_size):
        sub_aln = alignment[:, i:i+window_size]
        site_pis = []
        
        for col in range(window_size):
            column = sub_aln[:, col].upper()
            # Filter out gaps for pi calculation
            bases = [b for b in column if b in "ACGT"]
            if len(bases) < 2:
                site_pis.append(0)
                continue
            
            # Pi formula: 1 - sum(frequency^2)
            counts = {b: bases.count(b) for b in set(bases)}
            n = len(bases)
            pi_site = 1 - sum((count/n)**2 for count in counts.values())
            # Correction factor n/(n-1)
            site_pis.append(pi_site * (n / (n - 1)))
            
        pi_values.append(np.mean(site_pis))
    
    return pi_values

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--plot_data", required=True)
    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    
    # 1. Generate Consensus (threshold 0.5 means majority rule)
    consensus = summary.dumb_consensus(threshold=0.5)
    with open(args.consensus, "w") as f:
        f.write(f">Pangenome_Consensus\n{consensus}\n")

    # 2. Calculate Pi Diversity
    pi_data = calculate_pi(alignment)
    df = pd.DataFrame({"Position": range(len(pi_data)), "Pi": pi_data})
    df.to_csv(args.plot_data, sep="\t", index=False)

if __name__ == "__main__":
    main()
