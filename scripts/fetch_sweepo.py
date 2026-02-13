#!/usr/bin/env python3
import os
import csv
import argparse
import time
import sys
import socket
from Bio import Entrez, SeqIO

# --- Configuration for stability ---
socket.setdefaulttimeout(120)
Entrez.tool = "Vipan_Pangenomics_Final"

def search_nucleotide(query, email):
    Entrez.email = email
    print(f"[*] Querying NCBI...")
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=5000)
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        print(f"[!] Search Error: {e}")
        sys.exit(1)

def fetch_details(id_list, webenv, query_key, email, batch_size=30): # Smaller batch for stability
    Entrez.email = email
    total = len(id_list)
    records = []
    
    for i in range(0, total, batch_size):
        end = min(i + batch_size, total)
        print(f"  > Downloading records {i+1} to {end}...")
        attempt = 0
        while attempt < 3:
            try:
                handle = Entrez.efetch(
                    db="nucleotide", rettype="gb", retmode="text",
                    retstart=i, retmax=batch_size, webenv=webenv, query_key=query_key
                )
                batch_records = list(SeqIO.parse(handle, "genbank"))
                records.extend(batch_records)
                handle.close()
                break
            except Exception as e:
                attempt += 1
                print(f"  [!] Retry {attempt}/3: {e}")
                time.sleep(20)
        time.sleep(1)
    return records

def filter_and_save(records, min_len, max_len, fasta_out, tsv_out):
    final_records = []
    metadata = []
    
    os.makedirs(os.path.dirname(os.path.abspath(fasta_out)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(tsv_out)), exist_ok=True)

    for rec in records:
        desc = rec.description.lower()
        seq_len = len(rec.seq)
        # Avoid fragments but keep genomes
        is_partial = any(x in desc for x in ["partial", "fragment", "cds"])
        
        if not is_partial and (min_len <= seq_len <= max_len):
            final_records.append(rec)
            
            source = [f for f in rec.features if f.type == "source"]
            qual = source[0].qualifiers if source else {}
            
            # --- Country Extraction Logic ---
            country = qual.get("country", ["N/A"])[0]
            lat_lon = qual.get("lat_lon", ["N/A"])[0]
            isolate = qual.get("isolate", qual.get("strain", ["N/A"]))[0]

            # Logic: Infer Burkina Faso from GPS or Isolate name
            if country == "N/A":
                if lat_lon != "N/A" and ("N" in lat_lon or "W" in lat_lon):
                    country = "Check GPS (Possibly West Africa)"
                elif isolate.startswith("BF"):
                    country = "Burkina Faso (Inferred)"

            metadata.append({
                "accession": rec.id,
                "organism": rec.annotations.get("organism", "N/A"),
                "length": seq_len,
                "country": country,
                "lat_lon": lat_lon,
                "host": qual.get("host", ["N/A"])[0],
                "collection_date": qual.get("collection_date", ["N/A"])[0],
                "isolate": isolate
            })

    with open(fasta_out, "w") as f:
        SeqIO.write(final_records, f, "fasta")
    
    if metadata:
        keys = metadata[0].keys()
        with open(tsv_out, 'w', newline='', encoding='utf-8') as t:
            writer = csv.DictWriter(t, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            writer.writerows(metadata)

    return len(final_records)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", required=True)
    parser.add_argument("--min_len", type=int)
    parser.add_argument("--max_len", type=int)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--tsv", required=True)
    args = parser.parse_args()

    virus_taxa = ['"Sweet potato leaf curl virus"[Organism]', 
                  '"Ipomoea yellow vein virus"[All Fields]', 
                  '"Begomovirus ipomoeae"[All Fields]',
                  '"Begomovirus ipomoeageorgiaense"[All Fields]',
                  '"Begomovirus ipomoeaguangxiense"[All Fields]',
                  '"Begomovirus ipomoeahenanense"[All Fields]',
                  '"Begomovirus ipomoeahubeiense"[All Fields]',
                  '"Begomovirus ipomoeakoreaense"[All Fields]',
                  '"Begomovirus ipomoeamusivi"[All Fields]',
                  '"Begomovirus ipomoeasaopauloense"[All Fields]',
                  '"Begomovirus ipomoeashandongense"[All Fields]',
                  '"Begomovirus ipomoeasichuanprimi"[All Fields]',
                  '"Begomovirus ipomoeasichuansecundi"[All Fields]',
                  '"Begomovirus ipomoeasouthcarolinaense"[All Fields]',
                  '"Sweet potato leaf curl Lanzarote virus"[All Fields]', 
                  '"Sweet potato leaf curl Spain virus"[All Fields]']
    
    query = f"({' OR '.join(virus_taxa)}) AND {args.min_len}:{args.max_len}[Sequence Length] NOT (partial[Title] OR \"fragment\"[Title])"

    res = search_nucleotide(query, args.email)
    if int(res["Count"]) == 0:
        print("No records found."); open(args.fasta, 'w').close(); open(args.tsv, 'w').close()
        return

    raw = fetch_details(res["IdList"], res["WebEnv"], res["QueryKey"], args.email)
    kept = filter_and_save(raw, args.min_len, args.max_len, args.fasta, args.tsv)
    print(f"Done. Kept {kept} sequences.")

if __name__ == "__main__":
    main()
