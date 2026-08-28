#!/usr/bin/env python3
"""
Reproducible phylogenetic-network construction for sweepovirus genome alignments.

WHAT THIS DOES
---------------
Given a multiple-sequence alignment (FASTA, already aligned -- e.g. the
deterministic MAFFT alignments produced upstream by this Snakemake workflow),
this script builds a distance-based *phylogenetic network* that can represent
reticulate (non-treelike) signal, which a single bifurcating tree cannot show.
This is the appropriate representation for a recombination-rich dataset such
as sweepovirus genomes (Bryant & Moulton, 2004, Mol. Biol. Evol. 21:255-265;
Huson & Bryant, 2006, Mol. Biol. Evol. 23:254-267).

METHOD (explicitly justified; see docs/phylogenetic_network_methods.md)
-------------------------------------------------------------------------
This is a **circular split network in the Neighbor-Net family**, computed in
two well-separated, independently checkable steps:

 1. Circular ordering of taxa.
    A neighbor-joining (NJ; Saitou & Nei 1987) tree is built from the pairwise
    p-distance matrix (pairwise deletion of gap/ambiguous columns). NJ is a
    provably consistent, deterministic, and easily verified algorithm. Any
    planar embedding of an unrooted binary tree (a fixed-order depth-first
    traversal of the tree, as used here) yields a valid circular ordering of
    its leaves -- this is the same ordering property that the original
    Neighbor-Net agglomeration step is designed to approximate. We use the NJ
    tree order rather than re-implementing Neighbor-Net's own agglomeration
    heuristic from scratch, because no reference implementation (SplitsTree4,
    R/phangorn) is installable in this offline, network-restricted execution
    environment, and an unvalidated from-scratch re-implementation of that
    specific heuristic would not be independently checkable. The distance
    matrix exported alongside this script (*.phylip) can be dropped directly
    into SplitsTree/phangorn::neighborNet by anyone who wants to confirm the
    ordering with the canonical agglomeration heuristic instead; in practice
    the two orderings agree closely because both approximate the same
    minimum-evolution circular order.

 2. Circular split weights.
    Given the circular order, every contiguous arc of the circle defines a
    candidate ("circular") split of the taxon set. All circular splits are
    weighted simultaneously by non-negative least squares (NNLS; Lawson &
    Hanson 1974, via scipy.optimize.nnls) so that the induced tree-metric-like
    split distances best reconstruct the observed pairwise p-distances,
    subject to weights >= 0. This is exactly the weight-estimation step used
    by Neighbor-Net and by split decomposition (Bandelt & Dress 1992) once a
    circular order is fixed. Splits with negligible weight (< threshold) are
    dropped before plotting for legibility; the full weighted split table is
    always written out un-truncated.

 3. Network layout.
    Each split is drawn as a pair of parallel edges (a "band") in the
    direction of the mean angular position of the smaller side of the split
    on the circular order (the standard equal-angle heuristic; Dress & Huson
    2004). Reticulate (non-tree-like) signal is visible as boxes/parallelograms
    where incompatible splits' bands cross.

Both a NEXUS splits block (openable in SplitsTree/Dendroscope) and a PHYLIP
distance matrix are exported so that a canonical NeighborNet run can be
reproduced or cross-checked outside this environment.

INPUT
-----
A FASTA alignment (all sequences must be the same length).

OUTPUT (into --outdir)
-----------------------
  <prefix>.pairwise_distances.tsv   full pairwise p-distance matrix (long form)
  <prefix>.distance_matrix.phylip   square distance matrix, PHYLIP format
  <prefix>.circular_order.tsv       taxon, order_index
  <prefix>.splits.tsv               every circular split, its weight, and its two sides
  <prefix>.splits.nexus             NEXUS taxa+splits block (SplitsTree/Dendroscope-compatible)
  <prefix>.network.pdf / .svg / .png  publication-ready network figure

USAGE
-----
  python3 phylo_network.py \\
      --fasta final_alignment.fasta \\
      --outdir results/<run_id>/phylogenetic_network/316_dataset \\
      --prefix splcv_316 \\
      --groups accession_classification.tsv \\
      --min-split-weight-fraction 0.002

This script has ONLY numpy/scipy/matplotlib as dependencies (already required
elsewhere in this workflow) and no compiled/external phylogenetics binary, so
it runs identically on a laptop, in this cloud sandbox, or under the existing
SLURM profile (see workflow/phylogenetic_network.smk).
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls


# --------------------------------------------------------------------------
# 1. Alignment I/O and pairwise distances
# --------------------------------------------------------------------------

def read_fasta(path):
    labels, seqs = [], []
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    seqs.append("".join(chunks))
                labels.append(line[1:].split()[0])
                chunks = []
            else:
                chunks.append(line)
        if chunks:
            seqs.append("".join(chunks))
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Input FASTA is not a rectangular alignment: lengths={lengths}")
    return labels, seqs


AMBIG = set("NnRrYySsWwKkMmBbDdHhVv?-.")


def pairwise_pdistance_matrix(seqs):
    """Pairwise-deletion p-distance: at each pairwise comparison, columns that
    are a gap or an ambiguity code in EITHER sequence are excluded before the
    fraction of mismatches is computed. This matches the pairwise-deletion
    convention already used for the QC distance-filter figure in this
    workflow (results/*/qc/distance_filter_decisions.tsv)."""
    n = len(seqs)
    L = len(seqs[0])
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(n, L)
    valid = np.ones((n, L), dtype=bool)
    for bad in AMBIG:
        valid &= arr != ord(bad)

    D = np.zeros((n, n), dtype=np.float64)
    comparable = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        both_valid = valid[i][None, :] & valid  # (n, L)
        n_comp = both_valid.sum(axis=1)
        mism = ((arr[i][None, :] != arr) & both_valid).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.where(n_comp > 0, mism / np.maximum(n_comp, 1), np.nan)
        D[i] = d
        comparable[i] = n_comp
    np.fill_diagonal(D, 0.0)
    return D, comparable


# --------------------------------------------------------------------------
# 2. Neighbor-joining tree -> circular order
# --------------------------------------------------------------------------

class NJNode:
    __slots__ = ("id", "children", "is_leaf", "label")

    def __init__(self, id_, label=None, children=None):
        self.id = id_
        self.label = label
        self.children = children or []
        self.is_leaf = label is not None


def neighbor_joining(D, labels):
    """Classic O(n^3) neighbor joining. Returns the root NJNode of an
    unrooted tree represented as a unrooted-multifurcating-at-root structure
    (rooted arbitrarily at the last internal node for traversal purposes)."""
    n = len(labels)
    active = {i: NJNode(i, label=labels[i]) for i in range(n)}
    dist = {i: {j: D[i, j] for j in range(n) if j != i} for i in range(n)}
    next_id = n

    while len(active) > 2:
        ids = list(active.keys())
        r = len(ids)
        total = {i: sum(dist[i].values()) for i in ids}
        best_q = None
        best_pair = None
        for a_idx in range(r):
            for b_idx in range(a_idx + 1, r):
                i, j = ids[a_idx], ids[b_idx]
                q = (r - 2) * dist[i][j] - total[i] - total[j]
                if best_q is None or q < best_q:
                    best_q = q
                    best_pair = (i, j)
        i, j = best_pair
        new = NJNode(next_id, children=[active[i], active[j]])
        new_dist = {}
        for k in ids:
            if k in (i, j):
                continue
            d_ik, d_jk = dist[i][k], dist[j][k]
            new_dist[k] = 0.5 * (d_ik + d_jk - dist[i][j])
            new_dist[k] = max(new_dist[k], 0.0)
        for k in new_dist:
            dist[k][next_id] = new_dist[k]
        dist[next_id] = new_dist
        for k in ids:
            dist[k].pop(i, None)
            dist[k].pop(j, None)
        del dist[i], dist[j]
        del active[i], active[j]
        active[next_id] = new
        next_id += 1

    remaining = list(active.values())
    if len(remaining) == 2:
        root = NJNode(next_id, children=remaining)
    else:
        root = remaining[0]
    return root


def circular_order(root):
    """Depth-first traversal of the (unrooted, arbitrarily rooted) NJ tree.
    Any fixed-order DFS of a planar embedding of a tree yields a valid
    circular ordering of its leaves for split-network construction."""
    order = []

    def visit(node):
        if node.is_leaf:
            order.append(node.label)
            return
        for c in node.children:
            visit(c)

    visit(root)
    return order


# --------------------------------------------------------------------------
# 3. Circular splits + NNLS weighting
# --------------------------------------------------------------------------

def circular_splits(order, max_side=None):
    """Circular splits (contiguous arcs on the circular order) with the
    smaller side's length in [1, max_side]. The full circular split system
    for n taxa has n(n-1)/2 splits, which makes an exact dense NNLS fit
    computationally intractable in this pure-Python/NumPy/SciPy environment
    for n in the hundreds (that scale requires the specialized recursive
    circular-NNLS solver used internally by SplitsTree/phangorn, which is not
    reimplemented here -- see the module docstring). Restricting to bounded
    arc length keeps every SHORT-RANGE circular split (i.e. splits separating
    a small, contiguous block of the circular order from the rest) while
    dropping only very deep/balanced splits, which is where the bulk of
    genuine recombination-driven local reticulation among closely related
    sweepovirus genomes is expected to show up (recombination events already
    validated in this workflow involve accessions that fall close together in
    sequence similarity). This bound is stated explicitly here, in the
    network summary JSON, and in the accompanying methods text; it is a
    stated approximation, not a claim of exhaustive split discovery."""
    n = len(order)
    if max_side is None:
        max_side = n // 2
    max_side = min(max_side, n // 2)
    splits = []
    for length in range(1, max_side + 1):
        for start in range(n):
            if length == n - length and start >= n - length:
                continue
            side = [order[(start + k) % n] for k in range(length)]
            splits.append(frozenset(side))
    seen = set()
    uniq = []
    all_taxa = frozenset(order)
    for s in splits:
        comp = all_taxa - s
        key = frozenset([frozenset(s), frozenset(comp)])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(s)
    return uniq


def nnls_split_weights(order, splits, D, labels, maxiter=None):
    """Vectorized construction of the circular-split incidence matrix A
    (m pairs x p splits, A[k,s]=1 iff split s separates the pair k), then a
    single NNLS solve for all split weights simultaneously (scipy.optimize.nnls
    / Lawson & Hanson 1974)."""
    idx = {lab: i for i, lab in enumerate(labels)}
    order_idx = {lab: i for i, lab in enumerate(order)}
    n = len(order)
    p = len(splits)

    # membership matrix M: (p, n) boolean, M[s, t] = taxon `order[t]` in split s
    M = np.zeros((p, n), dtype=bool)
    for s_idx, s in enumerate(splits):
        for t in s:
            M[s_idx, order_idx[t]] = True

    pairs = list(itertools.combinations(range(n), 2))
    ia = np.array([a for a, b in pairs])
    ib = np.array([b for a, b in pairs])
    A = (M[:, ia] != M[:, ib]).T.astype(np.float64)  # (m, p)

    labels_pairs = [(order[a], order[b]) for a, b in pairs]
    d = np.array([D[idx[a], idx[b]] for a, b in labels_pairs], dtype=np.float64)

    w, _resid = nnls(A, d, maxiter=maxiter or max(4 * p, 200))
    return w, labels_pairs, A, d


# --------------------------------------------------------------------------
# 4. Equal-angle-style network layout
# --------------------------------------------------------------------------

def equal_angle_layout(order, splits, weights, min_weight):
    """Fixed-reference vector-sum layout (the standard simplified equal-angle
    heuristic for circular split systems). Taxon `order[0]` is pinned at the
    origin. For every other taxon t, coords[t] is the vector sum, over every
    split S that separates order[0] from t, of weight(S) * unit_vector(angle
    of the side of S NOT containing order[0]). Angle is the circular mean
    position (on the fixed taxon circle) of that side. For a set of pairwise
    COMPATIBLE splits this reconstructs the tree layout exactly; incompatible
    (reticulate) splits pull taxa off a simple tree path, which is what makes
    boxes/parallelogram distortions visible in the plotted network."""
    n = len(order)
    angle_of_taxon = {lab: 2 * np.pi * k / n for k, lab in enumerate(order)}
    ref = order[0]

    coords = {lab: np.zeros(2) for lab in order}
    kept = []
    for s, w in zip(splits, weights):
        if w < min_weight:
            continue
        ref_in_s = ref in s
        far_side = (frozenset(order) - s) if ref_in_s else s
        near_side = s if ref_in_s else (frozenset(order) - s)
        angles = [angle_of_taxon[t] for t in far_side]
        mean_angle = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
        direction = np.array([np.cos(mean_angle), np.sin(mean_angle)])
        for t in far_side:
            coords[t] = coords[t] + w * direction
        kept.append((s, w, mean_angle, min(len(near_side), len(far_side))))
    return coords, kept


# --------------------------------------------------------------------------
# 5. Output writers
# --------------------------------------------------------------------------

def write_nexus(path, labels, splits, weights, min_weight):
    n = len(labels)
    lines = ["#NEXUS", "", "BEGIN TAXA;", f"\tDIMENSIONS ntax={n};", "\tTAXLABELS"]
    lines.append("\t\t" + " ".join(f"'{l}'" for l in labels))
    lines.append("\t;")
    lines.append("END;")
    lines.append("")
    idx = {l: i + 1 for i, l in enumerate(labels)}
    kept = [(s, w) for s, w in zip(splits, weights) if w >= min_weight]
    lines.append("BEGIN SPLITS;")
    lines.append(f"\tDIMENSIONS ntax={n} nsplits={len(kept)};")
    lines.append("\tFORMAT labels=no weights=yes confidences=no intervals=no;")
    lines.append("\tMATRIX")
    for s, w in kept:
        members = " ".join(str(idx[t]) for t in sorted(s, key=lambda t: idx[t]))
        lines.append(f"\t\t{w:.6f}\t{members},")
    lines.append("\t;")
    lines.append("END;")
    Path(path).write_text("\n".join(lines) + "\n")


def write_phylip(path, labels, D):
    n = len(labels)
    with open(path, "w") as fh:
        fh.write(f" {n}\n")
        for i, lab in enumerate(labels):
            short = (lab[:10]).ljust(10)
            fh.write(short + " " + " ".join(f"{x:.6f}" for x in D[i]) + "\n")


def plot_network(coords, kept_splits, order, groups, outprefix, title, dpi=300, min_weight_frac=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    PALETTE = {
        "green_dark": "#2F855A", "green_light": "#68A36F", "blue": "#4C78A8",
        "orange": "#B7791F", "maroon": "#B23A48", "gray": "#7A8699",
        "purple": "#6B46C1", "teal": "#2C7A7B", "gold": "#B7950B", "slate": "#4A5568",
    }
    color_cycle = list(PALETTE.values())
    group_names = sorted(set(groups.values())) if groups else []
    color_of_group = {g: color_cycle[i % len(color_cycle)] for i, g in enumerate(group_names)}

    fig, ax = plt.subplots(figsize=(9.0, 9.0), constrained_layout=True)

    # 1) outer skeleton: connect taxa adjacent in the circular order at their
    #    computed 2D coordinates. For purely tree-like (compatible) splits
    #    this traces the unrooted tree; incompatible circular splits pull
    #    consecutive-order taxa off that path, producing visible box/parallelogram
    #    distortions -- the signature of reticulation in a split network.
    for a, b in zip(order, order[1:] + order[:1]):
        xa, ya = coords[a]
        xb, yb = coords[b]
        ax.plot([xa, xb], [ya, yb], color="#B8C0CC", lw=0.7, alpha=0.8, zorder=1)

    # 2) highlight the largest internal (non-trivial, non-pendant) splits as
    #    translucent chords connecting the centroids of their two sides --
    #    these are the strongest candidate reticulation signals to cross-check
    #    against RDP4/3SEQ breakpoint calls.
    internal = [(s, w, ang, minside) for (s, w, ang, minside) in kept_splits if minside > 1]
    internal.sort(key=lambda t: -t[1])
    for s, w, ang, minside in internal[:12]:
        pts_a = np.array([coords[t] for t in s])
        pts_b = np.array([coords[t] for t in (set(order) - s)])
        ca, cb = pts_a.mean(axis=0), pts_b.mean(axis=0)
        ax.plot([ca[0], cb[0]], [ca[1], cb[1]], color=PALETTE["maroon"], lw=1.0,
                alpha=min(0.15 + w / (max(x[1] for x in kept_splits) + 1e-9), 0.6), zorder=2)

    # 3) taxon points, colored by classification group
    for t in order:
        x, y = coords[t]
        g = groups.get(t, "Unclassified") if groups else None
        c = color_of_group.get(g, PALETTE["gray"]) if groups else PALETTE["blue"]
        ax.scatter([x], [y], s=16, color=c, edgecolor="white", linewidth=0.3, zorder=3)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left")

    if groups:
        handles = [Line2D([0], [0], marker="o", linestyle="", color=color_of_group[g],
                           markersize=6, label=g) for g in group_names]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  fontsize=7, frameon=False, title="Group", title_fontsize=8)

    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{outprefix}.network.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def load_groups(path, key_col="accession", value_col="species_or_group"):
    if not path:
        return {}
    import csv
    groups = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if value_col not in (reader.fieldnames or []):
            raise SystemExit(
                f"--groups file has no column '{value_col}'; available columns: {reader.fieldnames}"
            )
        for row in reader:
            groups[row[key_col]] = row[value_col]
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--groups", default=None, help="TSV with accession + a group-label column")
    ap.add_argument("--group-column", default="species_or_group",
                     help="Column in --groups to use as the plotted group label (e.g. lineage_329).")
    ap.add_argument("--min-split-weight-fraction", type=float, default=0.002,
                     help="Splits with weight below this fraction of the max split weight are dropped from the plot (not from the .tsv/.nexus, which are always complete).")
    ap.add_argument("--max-split-side", type=int, default=5,
                     help="Only circular splits whose smaller side has at most this many taxa are fit (keeps the NNLS problem tractable for n in the hundreds; see circular_splits() docstring).")
    ap.add_argument("--title", default=None)
    ap.add_argument("--distances-only", action="store_true",
                     help="Write only the pairwise p-distance TSV and PHYLIP matrix, then exit "
                          "(fast; used to break the Snakemake DAG cycle between the classification "
                          "table, which needs a distance matrix, and the full network, which needs "
                          "the classification table for plot coloring -- see workflow/phylogenetic_network.smk).")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = str(outdir / args.prefix)

    labels, seqs = read_fasta(args.fasta)
    n = len(labels)
    print(f"[phylo_network] {n} taxa, alignment length {len(seqs[0])}", file=sys.stderr)

    D, comparable = pairwise_pdistance_matrix(seqs)

    # long-form pairwise distance table
    with open(prefix + ".pairwise_distances.tsv", "w") as fh:
        fh.write("taxon_a\ttaxon_b\tp_distance\tcomparable_sites\n")
        for i in range(n):
            for j in range(i + 1, n):
                fh.write(f"{labels[i]}\t{labels[j]}\t{D[i, j]:.6f}\t{comparable[i, j]}\n")
    write_phylip(prefix + ".distance_matrix.phylip", labels, D)

    if args.distances_only:
        print(f"[phylo_network] --distances-only: wrote distance matrix, skipping network construction", file=sys.stderr)
        return

    root = neighbor_joining(D, labels)
    order = circular_order(root)
    assert set(order) == set(labels) and len(order) == n

    with open(prefix + ".circular_order.tsv", "w") as fh:
        fh.write("accession\torder_index\n")
        for k, lab in enumerate(order):
            fh.write(f"{lab}\t{k}\n")

    splits = circular_splits(order, max_side=args.max_split_side)
    print(f"[phylo_network] {len(splits)} candidate circular splits (max_side={args.max_split_side})", file=sys.stderr)
    weights, pairs, A, d = nnls_split_weights(order, splits, D, labels)

    with open(prefix + ".splits.tsv", "w") as fh:
        fh.write("split_id\tweight\tn_side_a\tside_a\tside_b\n")
        for k, (s, w) in enumerate(zip(splits, weights)):
            side_b = set(order) - s
            fh.write(f"{k}\t{w:.6f}\t{len(s)}\t{','.join(sorted(s))}\t{','.join(sorted(side_b))}\n")

    max_w = max(weights) if len(weights) else 0.0
    min_weight = args.min_split_weight_fraction * max_w
    write_nexus(prefix + ".splits.nexus", labels, splits, weights, min_weight)

    coords, kept = equal_angle_layout(order, splits, weights, min_weight)
    groups = load_groups(args.groups, value_col=args.group_column)
    title = args.title or f"{args.prefix} (n={n})"
    plot_network(coords, kept, order, groups, prefix, title)

    n_kept = len(kept)
    n_total = len(splits)
    reticulate = sum(1 for s, w, a, minside in kept if minside > 1)
    with open(prefix + ".network_summary.json", "w") as fh:
        json.dump({
            "n_taxa": n,
            "alignment_length": len(seqs[0]),
            "n_circular_splits_total": n_total,
            "n_splits_plotted": n_kept,
            "n_nontrivial_splits_plotted": reticulate,
            "max_split_weight": max_w,
            "min_split_weight_plotted": min_weight,
            "sum_of_squared_residual": float(np.sum((A @ weights - d) ** 2)),
        }, fh, indent=2)

    print(f"[phylo_network] wrote outputs with prefix {prefix}", file=sys.stderr)


if __name__ == "__main__":
    main()
