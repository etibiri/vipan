"""Synthetic, recombination-free backgrounds (coalescent tree, JC + gamma rates, 3000 sites)
with simulated recombinants at exact breakpoints. Same design as sim_make.py."""
import os
ROOT = os.environ.get("VIPAN_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")

import itertools, random
import numpy as np

OUT = f"{WORK}/data_syn"
N_REP, N_BG, NSITE = 30, 24, 3000
LENGTHS = [300, 800, 1500]
BINS = {"low": (0.03, 0.08), "mid": (0.08, 0.15), "high": (0.15, 0.30)}
TREE_HEIGHT = 0.16  # substitutions/site from root to tips
nprng = np.random.default_rng(20261007)
rng = random.Random(20261007)
B = np.array(list("ACGT"))

def coalescent(n):
    """Kingman tree -> list of (child, parent, branch_length); tips 0..n-1, ultrametric."""
    nodes = list(range(n)); nxt = n; t = 0.0; times = {i: 0.0 for i in nodes}
    edges = []
    while len(nodes) > 1:
        k = len(nodes)
        t += nprng.exponential(1.0 / (k * (k - 1) / 2))
        a, b = rng.sample(nodes, 2)
        nodes.remove(a); nodes.remove(b)
        edges += [(a, nxt, t - times[a]), (b, nxt, t - times[b])]
        times[nxt] = t; nodes.append(nxt); nxt += 1
    return edges, nxt - 1, t

def evolve(n, L):
    edges, root, height = coalescent(n)
    scale = TREE_HEIGHT / height
    rates = nprng.gamma(0.5, 2.0, L)
    seqs = {root: nprng.integers(0, 4, L)}
    for child, parent, bl in sorted(edges, key=lambda e: -e[1]):
        pass
    # parents have higher ids than children: process from root downwards
    for child, parent, bl in sorted(edges, key=lambda e: -e[1]):
        s = seqs[parent].copy()
        pchg = 0.75 * (1 - np.exp(-4.0 / 3.0 * bl * scale * rates))
        chg = nprng.random(L) < pchg
        s[chg] = (s[chg] + nprng.integers(1, 4, chg.sum())) % 4
        seqs[child] = s
    return ["".join(B[seqs[i]]) for i in range(n)]

def pdist(a, b):
    return sum(x != y for x, y in zip(a, b)) / len(a)

cells = list(itertools.product(BINS, LENGTHS))
truth = []
for r in range(N_REP):
    plan = [cells[(3 * r + k) % 9] for k in range(3)]
    for attempt in range(200):
        seqs = evolve(N_BG, NSITE)
        names = [f"S{r+1:02d}_{i+1:02d}" for i in range(N_BG)]
        A = np.array([list(s) for s in seqs])
        D = np.zeros((N_BG, N_BG))
        for i in range(N_BG):
            D[i] = (A != A[i]).mean(1)
        used, pairs, ok = set(), [], True
        for cls, L in plan:
            lo, hi = BINS[cls]
            cand = [(i, j) for i, j in itertools.permutations(range(N_BG), 2)
                    if i not in used and j not in used and lo <= D[i, j] < hi]
            if not cand:
                ok = False; break
            i, j = rng.choice(cand)
            s = rng.randrange(150, NSITE - L - 150)
            used.update([i, j]); pairs.append((cls, L, i, j, D[i, j], s))
        if ok:
            break
    else:
        raise SystemExit(f"rep {r+1}: no feasible tree")
    recs = []
    for k, (cls, L, i, j, d, s) in enumerate(pairs, 1):
        seq = list(seqs[i]); seq[s:s + L] = seqs[j][s:s + L]
        rid = f"SIM_r{r+1:02d}_e{k}"
        recs.append((rid, "".join(seq)))
        truth.append((r + 1, k, rid, names[i], names[j], s + 1, s + L, L, cls, round(d, 4)))
    with open(f"{OUT}/rep{r+1:02d}_ctrl.fa", "w") as f:
        for n_, s_ in zip(names, seqs):
            f.write(f">{n_}\n{s_}\n")
    with open(f"{OUT}/rep{r+1:02d}_test.fa", "w") as f:
        for n_, s_ in zip(names, seqs):
            f.write(f">{n_}\n{s_}\n")
        for rid, s_ in recs:
            f.write(f">{rid}\n{s_}\n")
with open(f"{OUT}/truth.tsv", "w") as f:
    f.write("rep\tevent\trec_id\tmajor\tminor\tstart\tend\tlength\tdiv_class\tparent_dist\n")
    for t in truth:
        f.write("\t".join(map(str, t)) + "\n")
print(len(truth), "events")
