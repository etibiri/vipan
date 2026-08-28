#!/usr/bin/env python3
"""
Cross-reference the circular split networks (329- and 316-sequence datasets)
against the classification table and the RDP4/3SEQ recombination results,
to produce the quantitative material needed for the manuscript's Materials
& Methods / Results / Discussion paragraphs on the phylogenetic-network
analysis.

Every number below is computed directly from the workflow's own output
files (no external or invented data):
  - results/.../phylogenetic_network/{329,316}_dataset/splcv_*.network_summary.json
  - results/.../phylogenetic_network/{329,316}_dataset/splcv_*.splits.tsv
  - results/.../phylogenetic_network/{329,316}_dataset/splcv_*.pairwise_distances.tsv
  - results/.../phylogenetic_network/accession_classification.tsv
"""
import json
import sys
from pathlib import Path

import math

import numpy as np
import pandas as pd

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")


def fisher_exact_greater(a, b, c, d):
    """Two-sided Fisher exact test p-value for a 2x2 table
    [[a, b], [c, d]], computed from scratch (no scipy dependency) via the
    hypergeometric distribution, summing all tables with the same margins
    whose probability is <= that of the observed table."""
    n = a + b + c + d
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d

    def log_choose(n_, k_):
        if k_ < 0 or k_ > n_:
            return float("-inf")
        return math.lgamma(n_ + 1) - math.lgamma(k_ + 1) - math.lgamma(n_ - k_ + 1)

    def log_prob(x):
        # P(a=x) under hypergeometric with margins row1,row2,col1,col2,n
        return (log_choose(col1, x) + log_choose(col2, row1 - x) - log_choose(n, row1))

    lo = max(0, row1 - col2)
    hi = min(row1, col1)
    obs_logp = log_prob(a)
    total = 0.0
    for x in range(lo, hi + 1):
        lp = log_prob(x)
        if lp <= obs_logp + 1e-9:
            total += math.exp(lp)
    odds = (a * d) / (b * c) if b * c != 0 else float("inf")
    return odds, min(1.0, total)


def load_dataset(tag):
    d = BASE / f"{tag}_dataset"
    prefix = f"splcv_{tag}"
    summary = json.load(open(d / f"{prefix}.network_summary.json"))
    splits = pd.read_csv(d / f"{prefix}.splits.tsv", sep="\t")
    dist = pd.read_csv(d / f"{prefix}.pairwise_distances.tsv", sep="\t")
    order = pd.read_csv(d / f"{prefix}.circular_order.tsv", sep="\t", header=None)
    return summary, splits, dist, order


def fit_statistic(summary, dist):
    """SplitsTree-style 'fit': 100 * (1 - SSR / sum(d_ij^2)), computed over
    the same set of pairs used by the NNLS fit (all n*(n-1)/2 pairs)."""
    ssq_d = float((dist["p_distance"] ** 2).sum())
    ssr = summary["sum_of_squared_residual"]
    return 100.0 * (1.0 - ssr / ssq_d), ssq_d


def taxon_side_sets(splits_df, order_labels):
    """Return list of (weight, set(side_a taxa), set(side_b taxa)) for
    non-trivial splits (both sides >1 taxon), sorted by weight desc."""
    all_taxa = set(order_labels)
    out = []
    for row in splits_df.itertuples(index=False):
        side_a = set(str(row.side_a).split(","))
        if row.n_side_a <= 1 or (len(all_taxa) - row.n_side_a) <= 1:
            continue
        side_b = all_taxa - side_a
        out.append((row.weight, side_a, side_b))
    out.sort(key=lambda t: -t[0])
    return out


def split_vs_lineage_conflict(side_a, side_b, lineage_of):
    """A split is 'clean' with respect to the lineage partition if each side
    is (near-)homogeneous in lineage membership; 'conflicting' if it mixes
    lineages roughly evenly on both sides (reticulation-compatible signal).
    Returns (purity_a, purity_b, n_lineages_split_across)."""
    def purity(side):
        labs = [lineage_of.get(t) for t in side if lineage_of.get(t) is not None]
        if not labs:
            return None, 0
        vc = pd.Series(labs).value_counts()
        return vc.iloc[0] / len(labs), len(vc)

    pa, na = purity(side_a)
    pb, nb = purity(side_b)
    labs_a = {lineage_of.get(t) for t in side_a if lineage_of.get(t) is not None}
    labs_b = {lineage_of.get(t) for t in side_b if lineage_of.get(t) is not None}
    shared = labs_a & labs_b
    return pa, pb, len(shared)


def main():
    cls = pd.read_csv(BASE / "accession_classification.tsv", sep="\t")
    cls_idx = cls.set_index("accession")

    report = {}

    for tag in ("329", "316"):
        summary, splits, dist, order = load_dataset(tag)
        order_labels = order[0].tolist() if order.shape[1] == 1 else order.iloc[:, 0].tolist()
        fit_pct, ssq_d = fit_statistic(summary, dist)

        lineage_col = f"lineage_{tag}"
        lineage_of = cls_idx[lineage_col].to_dict()

        nontrivial = taxon_side_sets(splits, order_labels)
        top20 = nontrivial[:20]

        # How many of the top-20 largest-weight non-trivial splits are
        # "clean" wrt the lineage partition (>=90% purity both sides) vs
        # "conflicting" (mix >=2 lineages substantially on both sides)?
        # NOTE: with n_lineages >= 15-21 in this dataset, the LARGER side of
        # a highly lopsided split (typical for these top-weight splits,
        # which isolate small clades) will essentially never reach >=90%
        # purity for a single lineage even for a perfectly tree-like split
        # -- that would require one lineage to hold >=90% of the entire
        # remaining taxon set. So purity of the smaller/minority side (the
        # biologically interpretable "clade being split off") is the
        # informative signal here, not a symmetric two-sided purity test.
        clean, conflicting, ambiguous = 0, 0, 0
        conflict_rows = []
        for w, sa, sb in top20:
            pa, pb, n_shared = split_vs_lineage_conflict(sa, sb, lineage_of)
            minority_p = pa if len(sa) <= len(sb) else pb
            if minority_p is None:
                ambiguous += 1
            elif minority_p >= 0.90:
                clean += 1
            else:
                conflicting += 1
                conflict_rows.append((w, round(pa, 2), round(pb, 2), len(sa), len(sb)))

        # Recombinant enrichment: are RDP4-robust recombinants found more
        # often on the minority/conflicting side of the largest non-trivial
        # splits than expected by chance? Build a 2x2 contingency table
        # across the top-20 splits' smaller side vs rest.
        recomb_set = set(cls_idx.index[cls_idx["rdp4_robust_recombinant"] == True])
        taxa_in_dataset = set(order_labels)
        n_recomb_total = len(recomb_set & taxa_in_dataset)

        small_side_union = set()
        for w, sa, sb in top20:
            smaller = sa if len(sa) <= len(sb) else sb
            small_side_union |= smaller
        n_small = len(small_side_union)
        n_small_recomb = len(small_side_union & recomb_set)
        n_rest = len(taxa_in_dataset) - n_small
        n_rest_recomb = n_recomb_total - n_small_recomb
        try:
            odds, pval = fisher_exact_greater(
                n_small_recomb, n_small - n_small_recomb,
                n_rest_recomb, n_rest - n_rest_recomb,
            )
        except Exception:
            odds, pval = float("nan"), float("nan")

        # Species / lineage composition summary
        in_ds = cls[cls[f"in_{tag}_dataset"] == True]
        species_counts = in_ds[f"species_or_group_{tag}"].value_counts().to_dict()
        lineage_counts = in_ds[lineage_col].value_counts().to_dict()
        n_lineages = in_ds[lineage_col].nunique()
        n_species_groups = in_ds[f"species_or_group_{tag}"].nunique()

        report[tag] = {
            "n_taxa": summary["n_taxa"],
            "alignment_length": summary["alignment_length"],
            "n_circular_splits_total": summary["n_circular_splits_total"],
            "n_splits_plotted": summary["n_splits_plotted"],
            "n_nontrivial_splits_plotted": summary["n_nontrivial_splits_plotted"],
            "fit_percent": round(fit_pct, 2),
            "sum_sq_residual": summary["sum_of_squared_residual"],
            "sum_sq_distances": round(ssq_d, 2),
            "n_species_or_groups": n_species_groups,
            "n_lineages": n_lineages,
            "species_counts_top10": dict(sorted(species_counts.items(), key=lambda x: -x[1])[:10]),
            "lineage_counts": lineage_counts,
            "top20_split_weight_range": [round(top20[-1][0], 5), round(top20[0][0], 5)] if top20 else None,
            "top20_clean_splits": clean,
            "top20_conflicting_splits": conflicting,
            "top20_ambiguous_splits": ambiguous,
            "conflicting_split_details": conflict_rows,
            "n_rdp4_robust_recombinants_in_dataset": n_recomb_total,
            "n_taxa_on_top20_minority_sides": n_small,
            "n_recombinants_on_top20_minority_sides": n_small_recomb,
            "fisher_odds_ratio": None if np.isnan(odds) else round(odds, 3),
            "fisher_pvalue": None if np.isnan(pval) else round(pval, 5),
        }

    print(json.dumps(report, indent=2, default=str))
    with open(BASE / "network_recombination_crossref.json", "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\n[wrote] {BASE / 'network_recombination_crossref.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
