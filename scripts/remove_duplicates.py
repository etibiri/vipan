#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import hashlib

def remove_duplicates(input_fasta, output_fasta):
    seen_hashes = set()
    kept = 0
    total = 0

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            # On calcule un hash unique basé sur la séquence nucléotidique
            # On met en majuscule pour éviter les faux doublons (a vs A)
            seq_hash = hashlib.md5(str(record.seq).upper().encode()).hexdigest()

            if seq_hash not in seen_hashes:
                seen_hashes.add(seq_hash)
                SeqIO.write(record, out_handle, "fasta")
                kept += 1
    
    return total, kept

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supprime les séquences 100% identiques")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    total, kept = remove_duplicates(args.input, args.output)
    print(f"Total traité : {total}")
    print(f"Séquences uniques gardées : {kept}")
    print(f"Doublons supprimés : {total - kept}")
