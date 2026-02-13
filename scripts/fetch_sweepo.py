#!/usr/bin/env python3
import os
import csv
import argparse
import time
import sys
import socket
from Bio import Entrez, SeqIO

# --- Configuration ---
# Increasing timeout for international NCBI connections
socket.setdefaulttimeout(90)
Entrez.tool = "Vipan_Pangenomics_Discovery"

def search_nucleotide(query, email):
    """Search NCBI and return the WebEnv/QueryKey for batch downloading."""
    Entrez.email = email
    print(f"[*] Querying NCBI Nucleotide with discovery mode...")
    try:
        # retmax=5000 to ensure we capture all potential virus accessions
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=5000)
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        print(f"[!] NCBI Search Error: {e}")
        sys.exit(1)

def fetch_details(id_list, webenv, query_key, email, batch_size=50):
    """Download GenBank records using batch processing."""
    Entrez.email = email
    total = len(id_list)
    records = []
    
    for i in range(0, total, batch_size):
        end = min(i + batch_size, total)
        print(f"  > Fetching records {i+1} to {end} of {total}...")
        attempt = 0
        while attempt < 3:
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=i,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key
                )
                batch_records = list(SeqIO.parse(handle, "genbank"))
                records.extend(batch_records)
                handle.close()
                break
            except Exception as e:
                attempt += 1
                print(f"  [!] Retry {attempt}/3 (Network lag?): {e}")
                time.sleep(15)
        time.sleep(1) # Respect NCBI servers
    return records

def filter_and_save(records, min_len, max_len, fasta_out, tsv_out):
    """Bioinformatics filtering based on length and exclusion of 'partial' tags."""
    final_records = []
    metadata = []
    rejected_count = 0

    # Create directories if they don't exist
    os.makedirs(os.path.dirname(os.path.abspath(fasta_out)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(tsv_out)), exist_ok=True)

    for rec in records:
        desc = rec.description.lower()
        seq_len = len(rec.seq)
        
        # WE RELAXED THE FILTER: 
        # We exclude only if 'partial' is explicitly mentioned.
        # We rely on the length (min_len to max_len) to define 'completeness'.
        is_partial = any(x in desc for x in ["partial", "fragment", "cds", "coat protein gene"])
        
        if not is_partial and (min_len <= seq_len <= max_len):
            final_records.append(rec)
            
            # Extracting source metadata
            source_features = [f for f in rec.features if f.type == "source"]
            qual = source_features[0].qualifiers if source_features else {}
            
            metadata.append({
                "accession": rec.id,
                "organism": rec.annotations.get("organism", "N/A"),
                "length": seq_len,
                "country": qual.get("country", ["N/A"])[0],
                "lat_lon": qual.get("lat_lon", ["N/A"])[0],
                "host": qual.get("host", ["N/A"])[0],
                "collection_date": qual.get("collection_date", ["N/A"])[0],
                "isolate": qual.get("isolate", qual.get("strain", ["N/A"]))[0]
            })
        else:
            rejected_count += 1

    # Save to FASTA
    with open(fasta_out, "w") as f_handle:
        SeqIO.write(final_records, f_handle, "fasta")
    
    # Save to TSV
    if metadata:
        keys = metadata[0].keys()
        with open(tsv_out, 'w', newline='', encoding='utf-8') as t_handle:
            dict_writer = csv.DictWriter(t_handle, fieldnames=keys, delimiter='\t')
            dict_writer.writeheader()
            dict_writer.writerows(metadata)

    return len(final_records), rejected_count

def main():
    parser = argparse.ArgumentParser(description="Vipan Viral Downloader - Discovery Mode")
    parser.add_argument("--email", required=True)
    parser.add_argument("--min_len", type=int, default=2000)
    parser.add_argument("--max_len", type=int, default=3500)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--tsv", required=True)
    
    args = parser.parse_args()

    # --- RELAXED DISCOVERY QUERY ---
    # We use size-based filtering directly in the NCBI query
    virus_taxa = [
        '"Sweet potato leaf curl virus"[Organism]',
        '"Ipomea yellow vein virus"[All Fields]',
        '"Ipomoea yellow vein virus"[All Fields]',
        '"Sweet potato leaf curl Lanzarote virus"[All Fields]',
        '"Sweet potato leaf curl Spain virus"[All Fields]'
    ]
    tax_query = " OR ".join(virus_taxa)
    
    # Query logic: (Taxonomy) AND (Length Range) NOT (Known partial markers)
    query = (
        f"({tax_query}) AND {args.min_len}:{args.max_len}[Sequence Length] "
        "NOT (partial[Title] OR \"partial sequence\"[Title] OR \"fragment\"[Title])"
    )

    search_res = search_nucleotide(query, args.email)
    id_list = search_res["IdList"]
    count = int(search_res["Count"])
    
    print(f"[+] NCBI Discovery found {count} records matching size and taxa.")

    if count == 0:
        open(args.fasta, 'a').close()
        open(args.tsv, 'a').close()
        return

    raw_records = fetch_details(id_list, search_res["WebEnv"], search_res["QueryKey"], args.email)
    kept, rejected = filter_and_save(raw_records, args.min_len, args.max_len, args.fasta, args.tsv)

    print(f"\n--- Results Summary ---")
    print(f"Total downloaded: {len(raw_records)}")
    print(f"Kept (in range):  {kept}")
    print(f"Rejected:         {rejected}")
    print(f"Files: {args.fasta}, {args.tsv}")

if __name__ == "__main__":
    main()
