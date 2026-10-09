#!/usr/bin/env python3
"""Project published breakpoints, given in the coordinates of a reference accession, onto the columns of the
final alignment.

Chain: GenBank position -> position in the oriented sequence (prepare_cassava.py may reverse-complement and
rotates every genome onto the origin nonanucleotide) -> column of aligned.fasta -> column of
alignment_final.fasta (postalign.py trims low-occupancy columns and writes the map as <final>.columns.tsv).

The rotation offset is not read from a log, it is recovered by matching the oriented sequence against the
raw one and its reverse complement, so a change in prepare_cassava.py cannot silently invalidate the
projection.
"""
import argparse

COMP = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path, wanted=None):
    out, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name and (wanted is None or name in wanted):
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.upper())
    if name and (wanted is None or name in wanted):
        out[name] = "".join(buf)
    return out


def rotation(raw, oriented):
    """Return (strand, offset) such that oriented == rot(seq, offset), seq being raw or its revcomp."""
    for strand, seq in (("+", raw), ("-", revcomp(raw))):
        if len(seq) != len(oriented):
            continue
        k = (seq + seq).find(oriented)
        if 0 <= k < len(seq):
            return strand, k
    raise SystemExit("oriented sequence is not a rotation of the raw sequence")


def raw_to_oriented(pos, strand, offset, L):
    """1-based position in the GenBank record -> 1-based position in the oriented sequence."""
    p = pos - 1
    if strand == "-":
        p = L - 1 - p
    return (p - offset) % L + 1


def ungapped_to_column(aligned_row):
    """1-based ungapped index -> 1-based alignment column."""
    m = {}
    k = 0
    for col, ch in enumerate(aligned_row, start=1):
        if ch != "-":
            k += 1
            m[k] = col
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--segment", required=True, choices=["A", "B"])
    ap.add_argument("--raw", required=True)
    ap.add_argument("--oriented", required=True)
    ap.add_argument("--aligned", required=True)
    ap.add_argument("--columns", required=True, help="<alignment_final>.columns.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.truth) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(head)}
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    rows = [r for r in rows if r[idx["segment"]] == a.segment]
    refs = {r[idx["ref_accession"]] for r in rows}

    raw = read_fasta(a.raw, refs)
    oriented = read_fasta(a.oriented, refs)
    aligned = read_fasta(a.aligned, refs)
    src_to_final = {}
    with open(a.columns) as fh:
        fh.readline()
        for line in fh:
            fin, src = line.split()
            src_to_final[int(src)] = int(fin)
    finals = sorted(src_to_final)

    cache = {}
    for ref in sorted(refs):
        if ref not in raw:
            continue
        strand, off = rotation(raw[ref], oriented[ref])
        cache[ref] = (strand, off, len(raw[ref]), ungapped_to_column(aligned[ref]))

    def project(ref, pos):
        strand, off, L, colmap = cache[ref]
        o = raw_to_oriented(int(pos), strand, off, L)
        src = colmap[o]
        if src in src_to_final:
            return src_to_final[src], "exact"
        nearest = min(finals, key=lambda c: abs(c - src))
        return src_to_final[nearest], f"nearest_kept(+{nearest - src})"

    with open(a.out, "w") as f:
        f.write("event_id\tsegment\tref_accession\tstrand\trotation_offset\tbegin_ref\tend_ref\t"
                "begin_col\tend_col\twraps_origin\tprojection\n")
        done = 0
        for r in rows:
            ref = r[idx["ref_accession"]]
            if ref not in cache:
                f.write("\t".join([r[idx["event_id"]], a.segment, ref, "NA", "NA", r[idx["begin_ref"]],
                                   r[idx["end_ref"]], "NA", "NA", r[idx["wraps_origin"]],
                                   "reference_absent"]) + "\n")
                continue
            b, qb = project(ref, r[idx["begin_ref"]])
            e, qe = project(ref, r[idx["end_ref"]])
            strand, off, L, _ = cache[ref]
            f.write("\t".join([r[idx["event_id"]], a.segment, ref, strand, str(off), r[idx["begin_ref"]],
                               r[idx["end_ref"]], str(b), str(e), r[idx["wraps_origin"]],
                               qb if qb == qe else f"{qb}/{qe}"]) + "\n")
            done += 1
    for ref, (strand, off, L, _) in cache.items():
        print(f"{ref}: strand {strand}, rotation offset {off}, length {L}")
    print(f"seg{a.segment}: {done}/{len(rows)} events projected -> {a.out}")


if __name__ == "__main__":
    main()
