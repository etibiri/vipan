#!/usr/bin/env python3
"""After MAFFT: trim columns carried by too few sequences, IUPAC codes -> N (OpenRDP, 3SEQ and GENECONV
reject them), drop sequences that are mostly gaps, collapse near-identical sequences (greedy), and cap the
number of sequences by farthest-point sampling.

Column trimming comes first and is not cosmetic: MAFFT opens insertion columns for a handful of divergent
records, so on the cassava DNA-A set every genome carried about 26% gaps over the full 3808-column width and
a 20% per-sequence gap filter emptied the alignment. Measured on the trimmed width, one sequence out of 974
is actually gappy.

Accessions listed in --keep are never removed, whatever the filter: they are the reference recombinants the
validation is scored against.

Writes the removed sequences to thinning_log.tsv and the trimmed-to-original column map to
<out>.columns.tsv, so breakpoints can be traced back to the untrimmed alignment."""
import argparse
import numpy as np


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0]); seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]


def read_keep(spec):
    if not spec:
        return set()
    try:
        with open(spec) as fh:
            toks = fh.read().replace(",", " ").split()
    except OSError:
        toks = spec.replace(",", " ").split()
    return {t.strip() for t in toks if t.strip() and not t.startswith("#")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--min-col-occ", type=float, default=0.1,
                    help="drop alignment columns present in fewer than this fraction of sequences")
    ap.add_argument("--max-gap-frac", type=float, default=0.2)
    ap.add_argument("--min-dist", type=float, default=0.005, help="drop a sequence closer than this to a kept one")
    ap.add_argument("--max-seqs", type=int, default=120, help="0 = no cap")
    ap.add_argument("--keep", default="", help="file or comma-separated list of accessions never to remove")
    a = ap.parse_args()
    names, seqs = read_fasta(a.aln)
    if not names:
        raise SystemExit("empty alignment")
    forced = read_keep(a.keep)
    missing = forced - set(names)
    arr = np.array([list(s) for s in seqs])
    arr[~np.isin(arr, list("ACGTN-"))] = "N"

    width0 = arr.shape[1]
    occ = (arr != "-").mean(0)
    col_keep = np.flatnonzero(occ >= a.min_col_occ)
    arr = arr[:, col_keep]
    ok = np.isin(arr, list("ACGT"))
    log = []

    is_forced = [n in forced for n in names]
    keep = []
    for i, n in enumerate(names):
        g = float((arr[i] == "-").mean())
        if g > a.max_gap_frac:
            if is_forced[i]:
                log.append((n, "gappy_but_forced", f"{g:.2f} gaps"))
                keep.append(i)
            else:
                log.append((n, "gappy", f"{g:.2f} gaps"))
        else:
            keep.append(i)

    def dist(i, j):
        m = ok[i] & ok[j]
        return float((arr[i][m] != arr[j][m]).mean()) if m.sum() > 500 else 1.0

    # forced accessions go first so they become the representatives of their near-identical clusters
    keep.sort(key=lambda i: (not is_forced[i], i))
    kept = []
    for i in keep:
        close = next((j for j in kept if dist(i, j) < a.min_dist), None)
        if close is None or is_forced[i]:
            if close is not None:
                log.append((names[i], "near_identical_but_forced", f"to {names[close]}"))
            kept.append(i)
        else:
            log.append((names[i], "near_identical", f"to {names[close]}"))
    if a.max_seqs and len(kept) > a.max_seqs:
        D = np.array([[dist(i, j) for j in kept] for i in kept])
        chosen = [k for k, i in enumerate(kept) if is_forced[i]]
        if not chosen:
            chosen = [int(D.sum(1).argmin())]
        if len(chosen) >= a.max_seqs:
            raise SystemExit(f"--keep lists {len(chosen)} sequences, more than --max-seqs={a.max_seqs}")
        mind = D[chosen].min(0)
        while len(chosen) < a.max_seqs:
            nxt = int(mind.argmax()); chosen.append(nxt); mind = np.minimum(mind, D[nxt])
        drop = set(range(len(kept))) - set(chosen)
        for k in sorted(drop):
            log.append((names[kept[k]], "capped", f"max_seqs={a.max_seqs}"))
        kept = [kept[k] for k in sorted(chosen)]
    kept.sort()
    with open(a.out, "w") as f:
        for i in kept:
            f.write(f">{names[i]}\n{''.join(arr[i])}\n")
    with open(a.out + ".columns.tsv", "w") as f:
        f.write("final_column\tsource_column\n")
        for k, c in enumerate(col_keep, start=1):
            f.write(f"{k}\t{c + 1}\n")
    with open(a.log, "w") as f:
        f.write("accession\tstatus\tdetail\n")
        f.write(f"-\tcolumns_trimmed\t{width0 - len(col_keep)} of {width0} below min_col_occ={a.min_col_occ}\n")
        for n in sorted(missing):
            f.write(f"{n}\tforced_absent\tlisted in --keep but not in the alignment\n")
        for r in log:
            f.write("\t".join(r) + "\n")
    print(f"postalign: {len(names)} in, {len(kept)} kept, {len(log)} removed; "
          f"width {width0} -> {len(col_keep)} columns"
          + (f"; {len(missing)} forced accessions absent" if missing else ""))


if __name__ == "__main__":
    main()
