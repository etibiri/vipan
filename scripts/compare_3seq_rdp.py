#!/usr/bin/env python3
"""
compare_3seq_rdp.py

Compare les evenements de recombinaison detectes par 3SEQ et par RDP5
(fichier "robust events min 3 methods") afin d'identifier les
recombinants confirmes par les deux approches.

Usage:
    python scripts/compare_3seq_rdp.py \
        --threeseq results/recombination/3seq_events_nonredundant.tsv \
        --rdp results/recombination/splcv_rdp_robust_events_min_3_methods.tsv \
        --out results/recombination/comparison_3seq_rdp.tsv \
        [--rdp-column recombinant_accessions]
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


def extract_accession(field: str) -> str:
    """
    Extrait le numero d'accession NCBI (ex: AB433786.1) en tete de
    chaine. Retourne une chaine vide si le champ est vide/absent.
    """
    if not field:
        return ""
    field = field.strip()
    if not field:
        return ""
    match = re.match(r"([A-Za-z]{1,4}_?\d+\.\d+)", field)
    if match:
        return match.group(1)
    parts = field.split()
    return parts[0] if parts else ""


def split_multi_accessions(field: str):
    """
    Certains fichiers RDP stockent plusieurs accessions dans un seul
    champ (ex: recombinant_accessions), separees par virgule,
    point-virgule ou barre verticale. Retourne la liste des
    accessions individuelles extraites.
    """
    if not field:
        return []
    field = field.strip()
    if not field:
        return []
    parts = re.split(r"[;,|]", field)
    accessions = [extract_accession(p) for p in parts]
    return [a for a in accessions if a]


def load_threeseq(path: str):
    """
    Retourne un dict accession_recombinant -> liste de lignes d'evenements 3SEQ.
    """
    events = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            recomb_acc = extract_accession(row["Recombinant"])
            if not recomb_acc:
                continue
            p1_acc = extract_accession(row["Best_parent_1"])
            p2_acc = extract_accession(row["Best_parent_2"])
            events[recomb_acc].append(
                {
                    "event_id": row["Event_ID"],
                    "parent_1": p1_acc,
                    "parent_2": p2_acc,
                    "breakpoint_1": row.get("Breakpoint_1", ""),
                    "breakpoint_2": row.get("Breakpoint_2", ""),
                    "p_value": row.get("P_value", ""),
                }
            )
    return events


def load_rdp(path: str, forced_column: str = None):
    """
    Retourne un dict accession_recombinant -> liste de lignes d'evenements RDP.
    Gere le cas ou une ligne contient plusieurs accessions recombinantes
    dans un seul champ (colonne 'recombinant_accessions').
    """
    events = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []

        if forced_column:
            recomb_col = forced_column
        else:
            # Priorite aux noms exacts les plus probables (evite de
            # matcher des colonnes de position comme recombinant_begin)
            preferred = ["recombinant_accessions", "recombinant_accession", "recombinant"]
            recomb_col = next((c for c in preferred if c in fieldnames), None)
            if recomb_col is None:
                recomb_col = next(
                    (c for c in fieldnames if "recombinant" in c.lower() and "accession" in c.lower()),
                    None,
                )
            if recomb_col is None:
                recomb_col = next((c for c in fieldnames if "recombinant" in c.lower()), fieldnames[0])

        print(f"[load_rdp] colonnes detectees : {fieldnames}", file=sys.stderr)
        print(f"[load_rdp] colonne utilisee comme recombinant : '{recomb_col}'", file=sys.stderr)

        skipped = 0
        multi_value_rows = 0
        for row in reader:
            accessions = split_multi_accessions(row.get(recomb_col, ""))
            if not accessions:
                skipped += 1
                continue
            if len(accessions) > 1:
                multi_value_rows += 1
            for acc in accessions:
                events[acc].append(row)

        if skipped:
            print(f"[load_rdp] {skipped} ligne(s) ignoree(s) (champ recombinant vide)", file=sys.stderr)
        if multi_value_rows:
            print(f"[load_rdp] {multi_value_rows} ligne(s) avec plusieurs accessions recombinantes", file=sys.stderr)
    return events, recomb_col


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threeseq", required=True, help="Fichier 3SEQ nonredundant events (TSV)")
    parser.add_argument("--rdp", required=True, help="Fichier RDP robust events min 3 methods (TSV)")
    parser.add_argument("--out", required=True, help="Fichier de sortie TSV")
    parser.add_argument(
        "--rdp-column",
        default=None,
        help="Force le nom de la colonne recombinant dans le fichier RDP (sinon auto-detection)",
    )
    args = parser.parse_args()

    threeseq_events = load_threeseq(args.threeseq)
    rdp_events, rdp_recomb_col = load_rdp(args.rdp, forced_column=args.rdp_column)

    threeseq_recombinants = set(threeseq_events.keys())
    rdp_recombinants = set(rdp_events.keys())

    shared = sorted(threeseq_recombinants & rdp_recombinants)
    only_3seq = sorted(threeseq_recombinants - rdp_recombinants)
    only_rdp = sorted(rdp_recombinants - threeseq_recombinants)

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "Recombinant_accession",
                "Detected_by",
                "N_3SEQ_events",
                "3SEQ_best_parent_1",
                "3SEQ_best_parent_2",
                "3SEQ_best_p_value",
                "N_RDP_events",
            ]
        )

        for acc in shared:
            best_3seq = min(threeseq_events[acc], key=lambda e: e["p_value"])
            writer.writerow(
                [
                    acc,
                    "3SEQ+RDP",
                    len(threeseq_events[acc]),
                    best_3seq["parent_1"],
                    best_3seq["parent_2"],
                    best_3seq["p_value"],
                    len(rdp_events[acc]),
                ]
            )

        for acc in only_3seq:
            best_3seq = min(threeseq_events[acc], key=lambda e: e["p_value"])
            writer.writerow(
                [
                    acc,
                    "3SEQ_only",
                    len(threeseq_events[acc]),
                    best_3seq["parent_1"],
                    best_3seq["parent_2"],
                    best_3seq["p_value"],
                    0,
                ]
            )

        for acc in only_rdp:
            writer.writerow([acc, "RDP_only", 0, "", "", "", len(rdp_events[acc])])

    total = len(threeseq_recombinants | rdp_recombinants)
    print(f"Colonne recombinant detectee dans le fichier RDP : '{rdp_recomb_col}'", file=sys.stderr)
    print(f"Recombinants 3SEQ            : {len(threeseq_recombinants)}", file=sys.stderr)
    print(f"Recombinants RDP (>=3 meth.) : {len(rdp_recombinants)}", file=sys.stderr)
    print(f"Confirmes par les deux       : {len(shared)}", file=sys.stderr)
    print(f"3SEQ uniquement              : {len(only_3seq)}", file=sys.stderr)
    print(f"RDP uniquement                : {len(only_rdp)}", file=sys.stderr)
    print(f"Total recombinants uniques    : {total}", file=sys.stderr)
    print(f"Resultats ecrits dans : {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
