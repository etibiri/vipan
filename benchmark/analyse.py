"""Score 3SEQ and GENECONV calls against the simulated truth.

Definitions (fixed before looking at results):
  flagged    : the simulated recombinant R is the child (3SEQ column C) or a member
               of a pair (GENECONV) with corrected p < 0.05.
  localised  : flagged, and the reported segment has Jaccard >= 0.3 with the true one.
  parents ok : among significant calls on R, the two reported parents are on opposite
               sides (one closer to the true major parent, one closer to the minor).
               3SEQ only; GENECONV does not name a recipient.
  partner ok : GENECONV partner of R is closer to the minor than to the major parent.
  consensus  : flagged and localised by 3SEQ and flagged and localised by GENECONV.
  new calls  : significant calls absent from the control run (same background, no
               simulated sequence) that involve background sequences only.
"""
import os
ROOT = os.environ.get("VIPAN_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")

import csv, glob, os, statistics as st
import numpy as np

SET = os.environ.get("SIM_SET", "")
O = f"{WORK}/out{SET}"
D = f"{WORK}/data{SET}"
ALPHA, JAC = 0.05, 0.3

# ---- distances between named sequences (jointly valid sites) ----
def read_fasta(p):
    names, seqs = [], []
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            names.append(l[1:].split()[0]); seqs.append([])
        elif l:
            seqs[-1].append(l.upper())
    return names, ["".join(s) for s in seqs]

def pdist(a, b):
    m = [(x, y) for x, y in zip(a, b) if x in "ACGT" and y in "ACGT"]
    if len(m) < 1000:
        return float("nan")
    return sum(x != y for x, y in m) / len(m)

def jac(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    inter = max(0, hi - lo + 1)
    union = (a[1] - a[0] + 1) + (b[1] - b[0] + 1) - inter
    return inter / union if union > 0 else 0.0

# ---- parsers ----
def parse_3seq(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        fh.readline()
        for line in fh:
            v = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(v) < 13:
                continue
            p = float(v[10])
            bps = []
            for loc in v[12:]:
                if "&" not in loc:
                    continue
                a, b = loc.split("&")
                a0 = [int(t) for t in a.strip().split("-")]
                b0 = [int(t) for t in b.strip().split("-")]
                bps.append((sum(a0) / len(a0), sum(b0) / len(b0)))
            if bps:
                rows.append(dict(P=v[0], Q=v[1], C=v[2], p=p, seg=bps[0]))
    return rows

def parse_gc(path):
    rows = []
    if not os.path.exists(path):
        return rows
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        if len(t) < 6 or ";" not in t[1]:
            continue
        a, b = t[1].split(";")[:2]
        if t[3] == ">":          # "> 1.0": corrected KA p-value not significant
            continue
        rows.append(dict(type=t[0], a=a, b=b, p=float(t[3]), seg=(int(t[4]), int(t[5]))))
    return rows

# ---- load truth ----
truth = list(csv.DictReader(open(f"{D}/truth.tsv"), delimiter="\t"))
by_rep = {}
for t in truth:
    t["start"], t["end"] = int(t["start"]), int(t["end"])
    by_rep.setdefault(int(t["rep"]), []).append(t)

res, fp_rows = [], []
for rep, evs in sorted(by_rep.items()):
    tag = f"rep{rep:02d}"
    names, seqs = read_fasta(f"{D}/{tag}_test.fa")
    S = dict(zip(names, seqs))
    sim_ids = {e["rec_id"] for e in evs}
    t3 = [r for r in parse_3seq(f"{O}/{tag}_test.3s.rec") if r["p"] < ALPHA]
    c3 = [r for r in parse_3seq(f"{O}/{tag}_ctrl.3s.rec") if r["p"] < ALPHA]
    tg = [r for r in parse_gc(f"{O}/{tag}_test.frags") if r["p"] < ALPHA]
    cg = [r for r in parse_gc(f"{O}/{tag}_ctrl.frags") if r["p"] < ALPHA]
    for e in evs:
        R, maj, mino = e["rec_id"], e["major"], e["minor"]
        truth_seg = (e["start"], e["end"])
        row = dict(rep=rep, rec=R, cls=e["div_class"], L=int(e["length"]),
                   d=float(e["parent_dist"]))
        # 3SEQ
        mine = [r for r in t3 if r["C"] == R]
        row["s3_flag"] = bool(mine)
        best = max(mine, key=lambda r: jac(r["seg"], truth_seg), default=None)
        row["s3_loc"] = bool(best and jac(best["seg"], truth_seg) >= JAC)
        if row["s3_loc"]:
            row["s3_err_start"] = abs(best["seg"][0] - e["start"])
            row["s3_err_end"] = abs(best["seg"][1] - e["end"])
        exact = relaxed = False
        for r in mine:
            if {r["P"], r["Q"]} == {maj, mino}:
                exact = True
            ds = []
            for x in (r["P"], r["Q"]):
                if x in S:
                    ds.append((pdist(S[x], S[maj]), pdist(S[x], S[mino])))
            if len(ds) == 2:
                a_maj = ds[0][0] < ds[0][1]; b_maj = ds[1][0] < ds[1][1]
                if a_maj != b_maj:
                    relaxed = True
        row["s3_par_exact"], row["s3_par_relaxed"] = exact, relaxed
        # GENECONV
        gm = [r for r in tg if R in (r["a"], r["b"])]
        row["gc_flag"] = bool(gm)
        gbest = max(gm, key=lambda r: jac(r["seg"], truth_seg), default=None)
        row["gc_loc"] = bool(gbest and jac(gbest["seg"], truth_seg) >= JAC)
        if row["gc_loc"]:
            row["gc_err_start"] = abs(gbest["seg"][0] - e["start"])
            row["gc_err_end"] = abs(gbest["seg"][1] - e["end"])
        pok = False
        for r in gm:
            x = r["b"] if r["a"] == R else r["a"]
            if x == mino:
                pok = True
            elif x in S and x not in sim_ids:
                if pdist(S[x], S[mino]) < pdist(S[x], S[maj]):
                    pok = True
        row["gc_partner_ok"] = pok
        row["cons"] = row["s3_loc"] and row["gc_loc"]
        res.append(row)
    # new calls on background only
    ctrl3_children = {r["C"] for r in c3}
    ctrl_pairs = {frozenset((r["a"], r["b"])) for r in cg}
    new3 = {r["C"] for r in t3
            if r["C"] not in sim_ids and r["P"] not in sim_ids and r["Q"] not in sim_ids
            and r["C"] not in ctrl3_children}
    newg = {frozenset((r["a"], r["b"])) for r in tg
            if r["a"] not in sim_ids and r["b"] not in sim_ids
            and frozenset((r["a"], r["b"])) not in ctrl_pairs}
    # consensus-new: 3SEQ child new and some new GENECONV pair containing it
    cons_new = {c for c in new3 if any(c in p for p in newg)}
    fp_rows.append(dict(rep=rep, ctrl3=len(ctrl3_children), ctrl_gc_pairs=len(ctrl_pairs),
                        new3=len(new3), new_gc=len(newg), new_cons=len(cons_new)))

# ---- summaries ----
def pct(x, n):
    return f"{100*x/n:.0f}% ({x}/{n})" if n else "NA"

def summarise(rows, label):
    n = len(rows)
    out = [label, f"  n events: {n}"]
    for k, nm in [("s3_flag", "3SEQ flagged"), ("s3_loc", "3SEQ localised"),
                  ("s3_par_exact", "3SEQ parents exact"), ("s3_par_relaxed", "3SEQ parents relaxed"),
                  ("gc_flag", "GENECONV flagged"), ("gc_loc", "GENECONV localised"),
                  ("gc_partner_ok", "GENECONV partner ok"), ("cons", "consensus localised")]:
        out.append(f"  {nm:22s} {pct(sum(r[k] for r in rows), n)}")
    for pre, nm in (("s3", "3SEQ"), ("gc", "GENECONV")):
        es = [r[f"{pre}_err_start"] for r in rows if f"{pre}_err_start" in r]
        ee = [r[f"{pre}_err_end"] for r in rows if f"{pre}_err_end" in r]
        if es:
            out.append(f"  {nm} breakpoint error (bp) median start {st.median(es):.0f} / end {st.median(ee):.0f}"
                       f"  (n={len(es)})")
    return "\n".join(out)

lines = [summarise(res, "ALL EVENTS")]
for cls in ("low", "mid", "high"):
    lines.append(summarise([r for r in res if r["cls"] == cls], f"DIVERGENCE {cls}"))
for L in (300, 800, 1500):
    lines.append(summarise([r for r in res if r["L"] == L], f"LENGTH {L}"))
n = len(fp_rows)
lines.append("CONTROL AND NEW CALLS (per replicate, 24 background genomes)")
lines.append(f"  control: median {st.median(r['ctrl3'] for r in fp_rows):.0f} genomes flagged as recombinant by 3SEQ; "
             f"median {st.median(r['ctrl_gc_pairs'] for r in fp_rows):.0f} GENECONV pairs")
for k, nm in [("new3", "3SEQ new background children"), ("new_gc", "GENECONV new background pairs"),
              ("new_cons", "consensus new")]:
    vals = [r[k] for r in fp_rows]
    lines.append(f"  {nm:32s} mean {sum(vals)/n:.2f}, replicates with >=1: {sum(v > 0 for v in vals)}/{n}")
text = "\n\n".join(lines)
print(text)
open(f"{O}/summary.txt", "w").write(text + "\n")

with open(f"{O}/per_event.tsv", "w") as f:
    keys = ["rep", "rec", "cls", "L", "d", "s3_flag", "s3_loc", "s3_par_exact", "s3_par_relaxed",
            "s3_err_start", "s3_err_end", "gc_flag", "gc_loc", "gc_partner_ok", "gc_err_start",
            "gc_err_end", "cons"]
    f.write("\t".join(keys) + "\n")
    for r in res:
        f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
with open(f"{O}/per_replicate_control.tsv", "w") as f:
    keys = list(fp_rows[0])
    f.write("\t".join(keys) + "\n")
    for r in fp_rows:
        f.write("\t".join(str(r[k]) for k in keys) + "\n")
