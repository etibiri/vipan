#!/usr/bin/env python3
"""Run native 3SEQ and GENECONV on an alignment and write one common event table.

events.tsv columns: tool, recombinant, parent_a, parent_b, start, end, p_corrected, note
  3SEQ     : recombinant = child (column C of the .3s.rec file); start/end = midpoints of the first listed
             breakpoint ranges (alignment columns, 1-based); p = Dunn-Sidak corrected.
  GENECONV : a significant shared fragment between two sequences; GENECONV does not say which one is the
             recipient, so recombinant = "NA"; p = Bonferroni-corrected Karlin-Altschul.
"""
import argparse, glob, os, subprocess, sys, time

ROOT = os.environ.get("VIPAN_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.environ.get("OPENRDP_BIN") or os.path.join(ROOT, "tools", "OpenRDP", "openrdp", "bin")


def run_3seq(aln, workdir):
    prefix = os.path.join(workdir, "3seq")
    for old in glob.glob(prefix + ".3s.*"):
        os.remove(old)  # 3SEQ asks before overwriting and then spins on an empty stdin
    t0 = time.time()
    subprocess.run([f"{BIN}/3Seq/3seq.Unix", "-full", aln, "-ptable", f"{BIN}/3Seq/myPvalueTable", "-d", "-id", prefix],
                   input=b"Y\n" * 20, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=workdir, timeout=6 * 3600)
    rows = []
    path = prefix + ".3s.rec"
    if os.path.exists(path):
        with open(path) as fh:
            fh.readline()
            for line in fh:
                v = [x.strip() for x in line.rstrip("\n").split("\t")]
                if len(v) < 13:
                    continue
                bps = []
                for loc in v[12:]:
                    if "&" in loc:
                        l, r = loc.split("&")
                        l = [int(t) for t in l.strip().split("-")]
                        r = [int(t) for t in r.strip().split("-")]
                        bps.append((sum(l) / len(l), sum(r) / len(r)))
                if bps:
                    rows.append(("3SEQ", v[2], v[0], v[1], round(bps[0][0]), round(bps[0][1]), float(v[10]), f"{len(bps)} breakpoint pairs"))
    return rows, time.time() - t0


def run_geneconv(aln, workdir):
    fa = os.path.join(workdir, "geneconv_in.fa")
    with open(aln) as src, open(fa, "w") as dst:
        dst.write(src.read())
    cfg = os.path.join(workdir, "geneconv.cfg")
    open(cfg, "w").write(f"#GCONV_CONFIG\n  -inputpath={fa}\n  -Indel_blocs\n  -Gscale=1\n  -Minlength=1\n"
                         "  -Minnpoly=2\n  -Minscore=2\n  -Maxoverlap=1\n")
    for old in (os.path.join(workdir, "geneconv_in.frags"), os.path.join(workdir, "geneconv_in.sum")):
        if os.path.exists(old):
            os.remove(old)
    t0 = time.time()
    subprocess.run([f"{BIN}/GENECONV/geneconv.Unix", f"-Seqfile={fa}", f"-Config={cfg}", "-nolog"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=workdir)
    rows = []
    path = os.path.join(workdir, "geneconv_in.frags")
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            if len(t) < 6 or ";" not in t[1] or t[3] == ">":
                continue
            a, b = t[1].split(";")[:2]
            p = float(t[3])
            if p < 0.05:
                rows.append(("GENECONV", "NA", a, b, int(t[4]), int(t[5]), p, t[0]))
    return rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--skip-geneconv", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    aln = os.path.abspath(a.aln)
    rows, t3 = run_3seq(aln, os.path.abspath(a.out))
    rows = [r for r in rows if r[6] < a.alpha]
    tg = 0.0
    if not a.skip_geneconv:
        g, tg = run_geneconv(aln, os.path.abspath(a.out))
        rows += g
    with open(f"{a.out}/events.tsv", "w") as f:
        f.write("tool\trecombinant\tparent_a\tparent_b\tstart\tend\tp_corrected\tnote\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    with open(f"{a.out}/timings.tsv", "w") as f:
        f.write(f"tool\tseconds\n3SEQ\t{t3:.1f}\nGENECONV\t{tg:.1f}\n")
    n3 = sum(r[0] == "3SEQ" for r in rows)
    print(f"3SEQ: {n3} significant triplets in {t3:.0f}s; GENECONV: {len(rows) - n3} fragments in {tg:.0f}s -> {a.out}/events.tsv")


if __name__ == "__main__":
    main()
