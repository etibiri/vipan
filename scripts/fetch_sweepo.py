import os
import csv
import argparse
import time
import sys
from Bio import Entrez, SeqIO

# --- Entrez Configuration ---
Entrez.tool = "Vipan_Snakemake_Downloader"

def search_nucleotide(query, email):
    """Search NCBI Nucleotide database and return search metadata."""
    Entrez.email = email
    print(f"[*] Querying NCBI: {query}")
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y")
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        print(f"[!] Critical Error during ESearch: {e}")
        sys.exit(1)

def fetch_details(id_list, webenv, query_key, email, batch_size=50):
    """Fetch GenBank records in batches using EFetch history."""
    Entrez.email = email
    total = len(id_list)
    records = []
    
    for i in range(0, total, batch_size):
        print(f"  > Downloading batch: {i}/{total}...")
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
                print(f"  [!] Retry {attempt}/3: {e}")
                time.sleep(5)
        time.sleep(0.5) # Compliance with NCBI rate limits
    return records

def filter_and_save(records, min_len, max_len, fasta_out, tsv_out):
    """Apply post-download filters and export data."""
    final_records = []
    metadata = []
    rejected_count = 0

    for rec in records:
        desc = rec.description.lower()
        seq_len = len(rec.seq)
        
        # Strict filtering criteria
        is_complete = "complete genome" in desc or "complete sequence" in desc
        is_partial = any(x in desc for x in ["partial", "gene", "cds", "coat protein"])
        
        if is_complete and not is_partial and (min_len <= seq_len <= max_len):
            final_records.append(rec)
            
            # Metadata extraction from 'source' feature
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
    with open(fasta_out, "w") as f:
        SeqIO.write(final_records, f, "fasta")
    
    # Save to TSV (Tab-separated)
    if metadata:
        keys = metadata[0].keys()
        with open(tsv_out, 'w', newline='', encoding='utf-8') as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
            dict_writer.writeheader()
            dict_writer.writerows(metadata)

    return len(final_records), rejected_count

def main():
    parser = argparse.ArgumentParser(description="SPLCV Genome Downloader for Snakemake")
    parser.add_argument("--email", required=True, help="User email for NCBI Entrez")
    parser.add_argument("--min_len", type=int, default=2500)
    parser.add_argument("--max_len", type=int, default=3000)
    parser.add_argument("--fasta", required=True, help="Output FASTA file path")
    parser.add_argument("--tsv", required=True, help="Output TSV metadata file path")
    
    args = parser.parse_args()

    # Default query for SPLCV complete genomes
    query = ('("Sweet potato leaf curl virus"[Organism] OR "sweet potato leaf curl virus"[All Fields]) '
             'AND (("complete genome"[Title] OR "complete genome"[All Fields]) OR '
             '("complete sequence"[Title] OR "complete sequence"[All Fields])) '
             'NOT (partial[Title] OR "partial sequence"[Title] OR "gene"[Title] OR '
             '"cds"[Title] OR "coat protein"[Title] OR "replication associated protein"[Title])')

    search_res = search_nucleotide(query, args.email)
    id_list = search_res["IdList"]
    count = int(search_res["Count"])
    
    if count == 0:
        print("[-] No records found.")
        # Create empty files to avoid Snakemake errors
        open(args.fasta, 'a').close()
        open(args.tsv, 'a').close()
        return

    raw_records = fetch_details(id_list, search_res["WebEnv"], search_res["QueryKey"], args.email)
    
    kept, rejected = filter_and_save(raw_records, args.min_len, args.max_len, args.fasta, args.tsv)

    print(f"\n--- Summary ---")
    print(f"Total NCBI hits: {count}")
    print(f"Validated: {kept} | Rejected: {rejected}")

if __name__ == "__main__":
    main()
