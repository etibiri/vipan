#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def rotate_and_flip(input_fasta, output_fasta, motif="TAATATTAC"):
    rotated_count = 0
    flipped_count = 0
    total = 0
    
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq_obj = record.seq.upper()
            seq_str = str(seq_obj)
            
            # 1. Try Forward Orientation
            pos = seq_str.find(motif)
            
            if pos != -1:
                new_seq = seq_str[pos:] + seq_str[:pos]
                rotated_count += 1
            else:
                # 2. Try Reverse Complement Orientation
                rc_seq_obj = seq_obj.reverse_complement()
                rc_seq_str = str(rc_seq_obj)
                pos_rc = rc_seq_str.find(motif)
                
                if pos_rc != -1:
                    new_seq = rc_seq_str[pos_rc:] + rc_seq_str[:pos_rc]
                    flipped_count += 1
                    rotated_count += 1
                else:
                    # 3. Motif not found in either orientation
                    new_seq = seq_str
            
            new_record = SeqRecord(
                Seq(new_seq),
                id=record.id,
                description=record.description
            )
            SeqIO.write(new_record, out_handle, "fasta")
    
    return total, rotated_count, flipped_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Normalize circular genomes: Flip RC and Rotate to motif")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--motif", default="TAATATTAC")
    args = parser.parse_args()

    total, rotated, flipped = rotate_and_flip(args.input, args.output, args.motif)
    print(f"Total processed: {total}")
    print(f"Normalized (Rotated): {rotated}")
    print(f"Flipped (Reverse Complement): {flipped}")
    print(f"Failed to find motif: {total - rotated}")
