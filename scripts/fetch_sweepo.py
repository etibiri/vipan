#!/usr/bin/env python3
import os
import csv
import argparse
import time
import sys
import socket
from Bio import Entrez, SeqIO

# --- Robustness Configuration ---
socket.setdefaulttimeout(120)  
Entrez.tool = "Vipan_Pangenomics_Final"

def search_nucleotide(query, email):
    """Search NCBI and handle the WebEnv session."""
    Entrez.email = email
    print(f"[*] Querying NCBI Nucleotide...")
    try:
        # Utilisation de post=True si la requête est très longue
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=5000)
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        print(f"[!] NCBI Search Error: {e}")
        sys.exit(1)

def fetch_details(id_list, webenv, query_key, email, batch_size=30):
    """Download records in small batches to avoid network drops."""
    Entrez.email = email
    total = len(id_list)
    records = []
    
    for i in range(0, total, batch_size):
        end = min(i + batch_size, total)
        print(f"  > Downloading batch {i//batch_size + 1}: {i+1} to {end}...")
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
                print(f"  [!] Retry {attempt}/3 (Network lag?): {e}")
                time.sleep(15)
        time.sleep(1) # Polite throttling
    return records

def infer_country(qualifiers, isolate):
    """
    Extracts country by checking /country, then /geo_loc_name, 
    then inferring from isolate or GPS.
    """
    # 1. Priorité au tag officiel /country
    country = qualifiers.get("country", [""])[0]
    
    # 2. Recherche dans /geo_loc_name (ex: "China: Henan")
    if not country or country == "N/A":
        geo_loc = qualifiers.get("geo_loc_name", [""])[0]
        if geo_loc:
            # On prend la partie avant les deux-points (ex: "China")
            country = geo_loc.split(":")[0].strip()

    # Si on a trouvé quelque chose de viable, on s'arrête
    if country and country != "N/A":
        # Nettoyage rapide (enlever les précisions après ":" si présentes dans /country)
        return country.split(":")[0].strip()

    # 3. Fallback: Inférence via le nom de l'isolat (codes pays courants)
    code_map = {
        "BR": "Brazil", "BF": "Burkina Faso", "US": "USA", "ES": "Spain", 
        "CN": "China", "KR": "South Korea", "TW": "Taiwan", "PE": "Peru"
    }
    for code, name in code_map.items():
        if f"[{code}:" in isolate or isolate.startswith(code) or f"-{code}" in isolate:
            return name

    # 4. Fallback: Analyse GPS rudimentaire
    lat_lon = qualifiers.get("lat_lon", [""])[0]
    if lat_lon:
        if "W" in lat_lon and "S" in lat_lon: return "South America (likely Brazil)"
        if "E" in lat_lon: return "Old World (Asia/Africa/Europe)"

    return "Unknown"

def filter_and_save(records, min_len, max_len, fasta_out, tsv_out):
    """Filter by length and exclude fragments, then extract metadata."""
    final_records = []
    metadata = []
    
    os.makedirs(os.path.dirname(os.path.abspath(fasta_out)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(tsv_out)), exist_ok=True)

    for rec in records:
        desc = rec.description.lower()
        seq_len = len(rec.seq)
        
        # Filtrage strict contre les séquences partielles
        is_partial = any(x in desc for x in ["partial", "fragment", "cds", "coat protein gene"])
        
        if not is_partial and (min_len <= seq_len <= max_len):
            final_records.append(rec)
            
            source_features = [f for f in rec.features if f.type == "source"]
            qual = source_features[0].qualifiers if source_features else {}
            
            isolate = qual.get("isolate", qual.get("strain", ["N/A"]))[0]
            
            # Utilisation de la nouvelle fonction d'inférence pays améliorée
            country = infer_country(qual, isolate)

            metadata.append({
                "accession": rec.id,
                "organism": rec.annotations.get("organism", "N/A"),
                "length": seq_len,
                "country": country,
                "lat_lon": qual.get("lat_lon", ["N/A"])[0],
                "host": qual.get("host", ["N/A"])[0],
                "collection_date": qual.get("collection_date", ["N/A"])[0],
                "isolate": isolate
            })

    # Sauvegarde FASTA
    with open(fasta_out, "w") as f_out:
        SeqIO.write(final_records, f_out, "fasta")
    
    # Sauvegarde TSV
    if metadata:
        with open(tsv_out, 'w', newline='', encoding='utf-8') as t_out:
            writer = csv.DictWriter(t_out, fieldnames=metadata[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(metadata)

    return len(final_records)

def main():
    parser = argparse.ArgumentParser(description="Vipan Viral Pangenomics Downloader")
    parser.add_argument("--email", required=True)
    parser.add_argument("--min_len", type=int, default=2200)
    parser.add_argument("--max_len", type=int, default=3300)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--tsv", required=True)
    args = parser.parse_args()

    virus_taxa = [
        '"Sweet potato leaf curl virus"[Organism]', 
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
        '"Sweet potato leaf curl Spain virus"[All Fields]',
        '"Begomovirus ipomoeacanaryense"[All Fields]',
        '"Begomovirus ipomoeachinaense"[All Fields]',
        '"sweet potato leaf curl Georgia virus"[All Fields]',
        '"sweet potato leaf curl Guangxi virus"[All Fields]',
        '"sweet potato leaf curl Henan virus"[All Fields]',
        '"sweet potato leaf curl Hubei virus"[All Fields]',
        '"sweet potato golden Korea vein virus"[All Fields]',
        '"sweet potato mosaic virus"[All Fields]',
        '"sweet potato leaf curl Sao Paulo virus"[All Fields]',
        '"sweet potato leaf curl Shandong virus"[All Fields]',
        '"sweet potato leaf curl Sichuan virus"[All Fields]',
        '"sweet potato leaf curl South Carolina virus"[All Fields]'
    ]
    
    # Construction de la requête
    query = (f"({' OR '.join(virus_taxa)}) "
             f"AND {args.min_len}:{args.max_len}[Sequence Length] "
             f"NOT (partial[Title] OR \"fragment\"[Title])")

    res = search_nucleotide(query, args.email)
    
    count = int(res["Count"])
    if count == 0:
        print("[-] No records found matching the query."); 
        open(args.fasta, 'w').close(); open(args.tsv, 'w').close()
        return

    print(f"[+] Found {count} accessions. Starting download...")
    raw = fetch_details(res["IdList"], res["WebEnv"], res["QueryKey"], args.email)
    kept = filter_and_save(raw, args.min_len, args.max_len, args.fasta, args.tsv)
    print(f"[*] Done. Kept {kept} sequences after filtering.")

if __name__ == "__main__":
    main()
