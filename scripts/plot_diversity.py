#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_pi(data_file, output_pdf):
    df = pd.read_csv(data_file, sep="\t")
    
    plt.figure(figsize=(12, 5))
    plt.plot(df['Position'], df['Pi'], color='#2c3e50', lw=1.5, label='Nucleotide Diversity (π)')
    plt.fill_between(df['Position'], df['Pi'], color='#3498db', alpha=0.3)
    
    # Generic Begomovirus Gene Annotations (Approximative positions)
    # Most Begomoviruses start at Ori (Pos 0)
    # CP is usually around 500-1100, Rep is around 1500-2600 (on reverse strand)
    genes = [
        (130, 480, 'V2', '#e74c3c'),
        (480, 1070, 'CP (V1)', '#f1c40f'),
        (1080, 1480, 'C3', '#2ecc71'),
        (1200, 1600, 'C2', '#9b59b6'),
        (1520, 2600, 'Rep (C1)', '#e67e22')
    ]
    
    for start, end, name, color in genes:
        plt.axvspan(start, end, alpha=0.15, color=color)
        plt.text((start + end) / 2, plt.gca().get_ylim()[1] * 0.9, name, 
                 horizontalalignment='center', fontweight='bold', color=color)

    plt.title('Pangenome Diversity Landscape (Sliding Window π)', fontsize=14)
    plt.xlabel('Genome Position (bp) - Starting at TAATATTAC', fontsize=12)
    plt.ylabel('Nucleotide Diversity (π)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim(0, df['Position'].max())
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_pdf)
    print(f"Plot saved to {output_pdf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    plot_pi(args.input, args.output)
