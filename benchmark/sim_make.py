"""Build simulated-recombinant datasets from the sweepovirus alignment.

Each replicate: 24 background genomes (control set) + 3 simulated recombinants
(test set). A recombinant is the major parent with one segment replaced by the
minor parent's sequence; coordinates are alignment columns (1-based, inclusive).
"""
import os
ROOT = os.environ.get("VIPAN_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")

import itertools, random, sys
import numpy as np

SRC = f"{WORK}/data/full316_sanitised.fa"
OUT = f"{WORK}/data"
N_REP, N_BG = 30, 24
LENGTHS = [300, 800, 1500]
BINS = {"low": (0.03, 0.08), "mid": (0.08, 0.15), "high": (0.15, 0.30)}

def read_fasta(p):
    names, seqs = [], []
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            names.append(l[1:].split()[0]); seqs.append([])
        elif l:
            seqs[-1].append(l.upper())
    return names, ["".join(s) for s in seqs]

names, seqs = read_fasta(SRC)
arr = np.array([list(s) for s in seqs])
valid = np.isin(arr, list("ACGT"))

def pdist(i, j):
    m = valid[i] & valid[j]
    if m.sum() < 2000:
        return None
    return float((arr[i][m] != arr[j][m]).mean())

cells = list(itertools.product(BINS, LENGTHS))  # 9 cells
rng = random.Random(20261006)
nseq, ncol = arr.shape
V = valid.astype(np.float32)
# pairwise p-distance over jointly valid sites (all 316 x 316, once)
codes = {c: (arr == c).astype(np.float32) for c in "ACGT"}
same = sum(codes[c] @ codes[c].T for c in "ACGT")
joint = V @ V.T
D = np.where(joint >= 2000, 1 - same / np.maximum(joint, 1), np.nan)
print("distance matrix done; median", np.nanmedian(D), flush=True)

def window_ok(a, b, L, tries=400):
    m = (valid[a] & valid[b]).astype(np.int32)
    cs = np.concatenate([[0], np.cumsum(m)])
    for _ in range(tries):
        s = rng.randrange(200, ncol - L - 200)
        if cs[s + L] - cs[s] >= 0.97 * L:
            return s
    return None

truth = []
for r in range(N_REP):
    plan = [cells[(3 * r + k) % 9] for k in range(3)]
    used, pairs = set(), []
    for cls, L in plan:
        lo, hi = BINS[cls]
        ia, ib = np.where((D >= lo) & (D < hi))
        idx = list(range(len(ia))); rng.shuffle(idx)
        found = None
        for t in idx[:3000]:
            a, b = int(ia[t]), int(ib[t])
            if a == b or a in used or b in used:
                continue
            s = window_ok(a, b, L)
            if s is not None:
                found = (a, b, float(D[a, b]), s); break
        if found is None:
            raise SystemExit(f"no pair for {cls} {L} rep {r+1}")
        a, b, d, s = found
        used.update([a, b]); pairs.append((cls, L, a, b, d, s))
    others = [i for i in range(nseq) if i not in used]
    bg = list(used) + rng.sample(others, N_BG - len(used))
    recs = []
    for k, (cls, L, a, b, d, s) in enumerate(pairs, 1):
        seq = list(seqs[a])
        seq[s:s + L] = seqs[b][s:s + L]
        rid = f"SIM_r{r+1:02d}_e{k}"
        recs.append((rid, "".join(seq)))
        truth.append((r + 1, k, rid, names[a], names[b], s + 1, s + L, L, cls, round(d, 4)))
    with open(f"{OUT}/rep{r+1:02d}_ctrl.fa", "w") as f:
        for i in bg:
            f.write(f">{names[i]}\n{seqs[i]}\n")
    with open(f"{OUT}/rep{r+1:02d}_test.fa", "w") as f:
        for i in bg:
            f.write(f">{names[i]}\n{seqs[i]}\n")
        for rid, s_ in recs:
            f.write(f">{rid}\n{s_}\n")

with open(f"{OUT}/truth.tsv", "w") as f:
    f.write("rep\tevent\trec_id\tmajor\tminor\tstart\tend\tlength\tdiv_class\tparent_dist\n")
    for t in truth:
        f.write("\t".join(map(str, t)) + "\n")
print(len(truth), "events,", N_REP, "replicates")
