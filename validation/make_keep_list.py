#!/usr/bin/env python3
"""Assign every cassava accession to a virus group, then build the list of accessions that postalign.py must
never remove: the two coordinate references of the truth table and a few representatives of every group that
appears in it, as recombinant or as parent. Without this the farthest-point cap can drop exactly the
sequences the validation is scored on.

Outputs, per segment: keep_seg<X>.txt (one accession per line) and groups_seg<X>.tsv (accession -> group).
"""
import argparse
import os
import re

# order matters: the first matching rule wins
RULES = [
    ("ACMBFV", lambda org, ctry, acc: "burkinafasoense" in org),
    ("CMMGV", lambda org, ctry, acc: "cassava mosaic madagascar" in org),
    ("EACMCV", lambda org, ctry, acc: "east african cassava mosaic cameroon" in org),
    ("EACMZV", lambda org, ctry, acc: "east african cassava mosaic zanzibar" in org),
    ("EACMMV", lambda org, ctry, acc: "east african cassava mosaic malawi" in org),
    ("EACMKV", lambda org, ctry, acc: "east african cassava mosaic kenya virus" in org),
    ("SACMV", lambda org, ctry, acc: "south african cassava mosaic" in org),
    ("EACMV-UG", lambda org, ctry, acc: org.startswith("east african cassava mosaic virus") and "uganda" in org),
    ("EACMV-KE", lambda org, ctry, acc: org.startswith("east african cassava mosaic virus")
     and ("-kenya" in org or "-ke2" in org)),
    ("EACMV", lambda org, ctry, acc: org.startswith("east african cassava mosaic virus")),
    ("ACMV", lambda org, ctry, acc: org == "begomovirus manihotis" or org.startswith("african cassava mosaic virus")),
    ("ICMV", lambda org, ctry, acc: "indian cassava mosaic" in org),
    ("SLCMV", lambda org, ctry, acc: "sri lankan cassava mosaic" in org),
]
EACMV_LIKE = {"EACMV", "EACMV-UG", "EACMV-KE", "EACMKV", "EACMZV", "EACMMV", "SACMV"}


def assign(org, ctry, acc):
    o, c = org.lower(), (ctry or "").lower()
    for name, rule in RULES:
        if rule(o, c, acc):
            return name
    return "other"


def read_manifest(path):
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(head)}
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) < len(head):
                continue
            rows.append({k: v[idx[k]] for k in idx})
    return rows


def fasta_names(path):
    return {l[1:].split()[0] for l in open(path) if l.startswith(">")}


def groups_in_truth(truth_path, segment):
    """Every group named in the truth table for this segment, as recombinant or as parent."""
    want, specific = set(), set()
    with open(truth_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(head)}
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if v[idx["segment"]] != segment or v[idx["scorable"]] != "yes":
                continue
            for col in ("recombinant_group", "minor_parent_group", "major_parent_group"):
                for g in v[idx[col]].split("|"):
                    g = g.strip()
                    if not g or g == "Unknown":
                        continue
                    if g == "EACMV-like":
                        want |= EACMV_LIKE
                    elif re.search(r"[\[\]:]", g):
                        specific.add(g)
                    else:
                        want.add(g)
            specific.add(v[idx["ref_accession"]])
    return want, specific


# isolate-level groups of the truth table that cannot be derived from the organism name
SPECIFIC_ACCESSIONS = {
    "EACMV-[SC]": ["JF909227.1", "JF909228.1"],          # the two Seychelles DNA-B sequences
    "EACMV-[AJ:Bam:AJ29AQ1:2009]": ["JF909211.1"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--oriented", required=True, help="oriented.fasta for this segment")
    ap.add_argument("--segment", required=True, choices=["A", "B"])
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--per-group", type=int, default=4, help="representatives forced per group")
    a = ap.parse_args()

    present = fasta_names(a.oriented)
    rows = read_manifest(a.manifest)
    gmap = {}
    for r in rows:
        if r["accession"] in present:
            gmap[r["accession"]] = assign(r["organism"], r.get("country", ""), r["accession"])

    want, specific = groups_in_truth(a.truth, a.segment)
    keep = []
    for g in sorted(want):
        members = sorted(acc for acc, gg in gmap.items() if gg == g)
        keep.extend(members[: a.per_group])
    for s in sorted(specific):
        for acc in SPECIFIC_ACCESSIONS.get(s, [s]):
            if acc in present:
                keep.append(acc)
    keep = sorted(set(keep))

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, f"keep_seg{a.segment}.txt"), "w") as f:
        f.write("\n".join(keep) + "\n")
    with open(os.path.join(a.out_dir, f"groups_seg{a.segment}.tsv"), "w") as f:
        f.write("accession\tgroup\n")
        for acc in sorted(gmap):
            f.write(f"{acc}\t{gmap[acc]}\n")

    counts = {}
    for g in gmap.values():
        counts[g] = counts.get(g, 0) + 1
    missing = [g for g in sorted(want) if counts.get(g, 0) == 0]
    print(f"seg{a.segment}: {len(gmap)} sequences, groups " +
          ", ".join(f"{g}={n}" for g, n in sorted(counts.items(), key=lambda x: -x[1])))
    print(f"seg{a.segment}: {len(keep)} accessions forced" +
          (f"; groups absent from the dataset: {', '.join(missing)}" if missing else ""))


if __name__ == "__main__":
    main()
