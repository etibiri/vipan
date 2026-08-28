#!/usr/bin/env python3
"""
Build the accession -> species / phylogenetic-cluster / genetic-group /
recombinant-status classification table used to annotate the phylogenetic
networks (scripts/build_phylogenetic_network.py) and cross-reference them
against the RDP4/3SEQ recombination results already produced by this
workflow.

CLASSIFICATION LOGIC (every step is deterministic and documented; nothing is
invented -- see docs/phylogenetic_network_methods.md for the full write-up)
---------------------------------------------------------------------------
1. species: the GenBank/NCBI "organism" field for each accession (data/raw
   metadata TSV), after two purely mechanical normalizations: (a) stripping a
   trailing "-[country:isolate:year]"-style bracket annotation, and (b)
   normalizing the "Sweetpotato"/"Sweet potato" spacing variant used
   inconsistently by different GenBank submitters for the same species name.
   No species names are merged or reassigned beyond this.

2. coarse_cluster: average-linkage hierarchical clustering of the pairwise
   p-distance matrix (independently, per dataset), cut at
   --coarse-threshold (default 0.09, i.e. ~91% identity -- the demarcation
   value commonly cited for begomovirus species; Brown et al., 2015,
   Archives of Virology 160:1593-1619; Muhire et al., 2013 SDT). This is a
   purely data-driven partition, cross-checked against the submitted species
   names (below) as a validation step, not the other way around.

3. species_or_group: for each coarse cluster, the most frequent normalized
   species name among its members, PROVIDED it covers >=70% of the cluster
   ("purity"); otherwise the cluster is reported as "Mixed cluster <id>"
   rather than being forced under a single (possibly wrong) species label.

4. lineage: coarse clusters larger than --split-size x (default 25x the
   split unit, i.e. >25 members) are further partitioned into internally
   coherent, roughly balanced sub-groups by average-linkage clustering
   restricted to that cluster's own distance submatrix (ceil(n/25)
   sub-groups). This exists specifically because the dominant
   "Sweet potato leaf curl virus" GenBank label is known to be
   non-monophyletic in this dataset (consistent with Mauriello et al. 2026);
   sub-groups are labelled "<species> genetic group <k>" -- never
   "subspecies", which is not a recognized rank for this genus.

5. Recombination status: rdp4_robust_recombinant / n_rdp4_robust_events_as_
   recombinant come directly from this workflow's own
   recombination/rdp_events_robust.tsv (the "recombinant_accessions" column,
   already filtered to statistically robust RDP4 events upstream).
   strict_validation_tier / strict_validation_class come directly from
   validation/event_validation_calls_strict.tsv (Tier A/B phylogenetic
   validation already computed by this workflow). No new recombination
   detection is performed here -- this script only cross-references
   already-validated calls against the network/clustering output.
"""
import argparse
import math
import re

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def normalize_species(raw):
    s = str(raw).strip()
    s = re.sub(r"\s*-?\s*\[[^\]]*\]\s*$", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"^Sweetpotato\b", "Sweet potato", s)
    return s


def load_distance_matrix(path):
    df = pd.read_csv(path, sep="\t")
    taxa = sorted(set(df.taxon_a) | set(df.taxon_b))
    idx = {t: i for i, t in enumerate(taxa)}
    n = len(taxa)
    D = np.zeros((n, n))
    for row in df.itertuples(index=False):
        i, j = idx[row.taxon_a], idx[row.taxon_b]
        D[i, j] = D[j, i] = row.p_distance
    return taxa, idx, D


def cluster_dataset(tag, dist_path, meta, coarse_t, split_size, min_to_split=15):
    taxa, idx, D = load_distance_matrix(dist_path)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")
    coarse = fcluster(Z, t=coarse_t, criterion="distance")
    coarse_df = pd.DataFrame({"accession": taxa, "coarse_cluster": coarse})

    m = coarse_df.merge(meta[["accession", "species"]], on="accession", how="left")
    group_label = {}
    for cl, sub in m.groupby("coarse_cluster"):
        vc = sub["species"].value_counts()
        top_species, top_n = vc.index[0], vc.iloc[0]
        purity = top_n / len(sub)
        if top_species == "unidentified unclassified sequences." and len(vc) > 1:
            top_species, top_n = vc.index[1], vc.iloc[1]
            purity = top_n / len(sub)
        group_label[cl] = top_species if purity >= 0.70 else f"Mixed cluster {int(cl)}"
    m["species_or_group"] = m["coarse_cluster"].map(group_label)

    lineage = {}
    for cl, sub in m.groupby("coarse_cluster"):
        accs = sub["accession"].tolist()
        if len(accs) <= min_to_split:
            for a in accs:
                lineage[a] = group_label[cl]
            continue
        sub_idx = [idx[a] for a in accs]
        sub_D = D[np.ix_(sub_idx, sub_idx)]
        sub_condensed = squareform(sub_D, checks=False)
        sub_Z = linkage(sub_condensed, method="average")
        k = max(2, math.ceil(len(accs) / split_size))
        sub_cl = fcluster(sub_Z, t=k, criterion="maxclust")
        base_name = group_label[cl]
        for a, sc in zip(accs, sub_cl):
            lineage[a] = f"{base_name} genetic group {sc}"
    m["lineage"] = m["accession"].map(lineage)
    return m.rename(columns={
        "coarse_cluster": f"coarse_cluster_{tag}",
        "species_or_group": f"species_or_group_{tag}",
        "lineage": f"lineage_{tag}",
    })[["accession", f"coarse_cluster_{tag}", f"species_or_group_{tag}", f"lineage_{tag}"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, help="data/raw/*.tsv with accession, organism columns")
    ap.add_argument("--distances-329", required=True)
    ap.add_argument("--distances-316", required=True)
    ap.add_argument("--rdp-robust", required=True)
    ap.add_argument("--strict-calls", required=True)
    ap.add_argument("--coarse-threshold", type=float, default=0.09)
    ap.add_argument("--split-size", type=int, default=25)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    meta = pd.read_csv(args.metadata, sep="\t")
    meta["species"] = meta["organism"].apply(normalize_species)

    r329 = cluster_dataset("329", args.distances_329, meta, args.coarse_threshold, args.split_size)
    r316 = cluster_dataset("316", args.distances_316, meta, args.coarse_threshold, args.split_size)

    master = meta.merge(r329, on="accession", how="left").merge(r316, on="accession", how="left")
    master["in_329_dataset"] = master["coarse_cluster_329"].notna()
    master["in_316_dataset"] = master["coarse_cluster_316"].notna()

    rdp = pd.read_csv(args.rdp_robust, sep="\t")
    recomb_counts = {}
    for _, row in rdp.iterrows():
        for a in str(row["recombinant_accessions"]).split(";"):
            a = a.strip()
            if a and a != "nan":
                recomb_counts[a] = recomb_counts.get(a, 0) + 1
    master["n_rdp4_robust_events_as_recombinant"] = master["accession"].map(recomb_counts).fillna(0).astype(int)
    master["rdp4_robust_recombinant"] = master["n_rdp4_robust_events_as_recombinant"] > 0

    strict = pd.read_csv(args.strict_calls, sep="\t")
    tier_map = strict.groupby("Recombinant_accession")["Evidence_tier"].apply(lambda x: ",".join(sorted(set(x)))).to_dict()
    class_map = strict.groupby("Recombinant_accession")["Strict_validation_class"].apply(lambda x: ",".join(sorted(set(x)))).to_dict()
    master["strict_validation_tier"] = master["accession"].map(tier_map)
    master["strict_validation_class"] = master["accession"].map(class_map)

    cols = ["accession", "organism", "species",
            "in_329_dataset", "in_316_dataset",
            "species_or_group_329", "lineage_329", "coarse_cluster_329",
            "species_or_group_316", "lineage_316", "coarse_cluster_316",
            "n_rdp4_robust_events_as_recombinant", "rdp4_robust_recombinant",
            "strict_validation_tier", "strict_validation_class",
            "host", "country", "collection_date", "isolate", "length"]
    master[cols].to_csv(args.output, sep="\t", index=False)
    print(f"[classify_accessions] wrote {args.output} ({len(master)} accessions)")


if __name__ == "__main__":
    main()
