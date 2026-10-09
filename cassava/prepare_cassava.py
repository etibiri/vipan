#!/usr/bin/env python3
"""Select one segment, give every genome the same start (circular origin motif TAATATTAC, either strand),
and drop exact duplicates. Writes oriented.fasta and prepare_log.tsv."""
import argparse, csv, os, re

MOTIF = "TAATATTAC"
COMP = str.maketrans("ACGTacgtRYKMBDHVrykmbdhv", "TGCAtgcaYRMKVHDByrmkvhdb")


def read_fasta(path):
    rec, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            rec[name] = []
        elif name and line:
            rec[name].append(line.upper().replace("U", "T"))
    return {k: "".join(v) for k, v in rec.items()}


def revcomp(s):
    return s.translate(COMP)[::-1]


def motif_hit(seq, max_mm=1):
    """Best circular placement of the motif on this strand: (mismatches, start) or None."""
    ext = seq + seq[:len(MOTIF) - 1]
    best = None
    for i in range(len(seq)):
        mm = sum(x != y for x, y in zip(ext[i:i + len(MOTIF)], MOTIF))
        if mm <= max_mm and (best is None or mm < best[0]):
            best = (mm, i)
            if mm == 0:
                break
    return best


def kmers(seq, k=12):
    ext = seq + seq[:k - 1]
    return {ext[i:i + k] for i in range(len(seq))}


def strand_calls(seqs, k=12):
    """Decide the strand of every record by shared k-mers with the other records, not by the motif:
    the origin stem-loop is nearly palindromic, so the motif can occur on both strands of one genome.
    Pass 1 uses one reference; pass 2 uses all confidently oriented records (each excluded from its own score)."""
    names = list(seqs)
    ref = None
    for n in names:
        f, r = motif_hit(seqs[n], 0), motif_hit(revcomp(seqs[n]), 0)
        if (f is None) != (r is None):
            ref = (n, seqs[n] if f else revcomp(seqs[n]))
            break
    if ref is None:
        n = names[0]
        ref = (n, seqs[n])
    from collections import Counter
    pools = Counter(kmers(ref[1], k))
    own = {n: (kmers(seqs[n], k), kmers(revcomp(seqs[n]), k)) for n in names}

    def score(n, pool, exclude):
        res = []
        for ks in own[n]:
            res.append(sum(1 for x in ks if pool.get(x, 0) - (1 if exclude and x in exclude else 0) > 0))
        return res

    calls = {}
    for n in names:
        sf, sr = score(n, pools, None)
        calls[n] = ("+" if sf >= sr else "-", sf, sr)
    for _ in range(2):  # pass 2: pool built from confident calls, then re-evaluate
        pool = Counter()
        member = {}
        for n in names:
            st, sf, sr = calls[n]
            if max(sf, sr) >= 5 * (min(sf, sr) + 1):
                ks = own[n][0] if st == "+" else own[n][1]
                pool.update(ks)
                member[n] = ks
        newcalls = {}
        for n in names:
            sf, sr = [sum(1 for x in ks if pool.get(x, 0) - (1 if n in member and x in member[n] else 0) > 0) for ks in own[n]]
            newcalls[n] = ("+" if sf >= sr else "-", sf, sr)
        calls = newcalls
    return calls


def orient(seq, strand):
    """Orient with the given strand, then rotate to the best motif placement (<=1 mismatch)."""
    s = seq if strand == "+" else revcomp(seq)
    hit = motif_hit(s)
    if hit is None:
        return None, "motif_absent"
    return s[hit[1]:] + s[:hit[1]], f"strand={strand};motif_mismatches={hit[0]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--segment", default="A", choices=["A", "B", "unassigned", "all"])
    ap.add_argument("--min-len", type=int, default=2500)
    ap.add_argument("--max-len", type=int, default=2900)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    seqs = read_fasta(a.raw)
    meta = {r["accession"]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}
    log, kept, seen = [], {}, {}
    cand = {}
    for acc, s in seqs.items():
        seg = meta.get(acc, {}).get("segment", "unassigned")
        if (a.segment == "all" or seg == a.segment) and a.min_len <= len(s) <= a.max_len:
            cand[acc] = s
    calls = strand_calls(cand) if cand else {}
    for acc, s in seqs.items():
        m = meta.get(acc, {})
        seg = m.get("segment", "unassigned")
        if a.segment != "all" and seg != a.segment:
            log.append((acc, "skipped", f"segment={seg}")); continue
        if not (a.min_len <= len(s) <= a.max_len):
            log.append((acc, "skipped", f"length={len(s)}")); continue
        st, sf, sr = calls[acc]
        if max(sf, sr) < 3 * (min(sf, sr) + 1):
            log.append((acc, "dropped", f"strand_uncertain(fwd={sf},rev={sr})")); continue
        o, why = orient(s, st)
        if o is None:
            log.append((acc, "dropped", why)); continue
        if o in seen:
            log.append((acc, "duplicate", f"identical to {seen[o]}")); continue
        seen[o] = acc
        kept[acc] = o
        log.append((acc, "kept", f"strand={why}"))
    with open(f"{a.out}/oriented.fasta", "w") as f:
        for acc, s in kept.items():
            f.write(f">{acc}\n{s}\n")
    with open(f"{a.out}/prepare_log.tsv", "w") as f:
        f.write("accession\tstatus\tdetail\n")
        for r in log:
            f.write("\t".join(r) + "\n")
    from collections import Counter
    print("prepare:", dict(Counter(r[1] for r in log)), "->", f"{a.out}/oriented.fasta")


if __name__ == "__main__":
    main()
