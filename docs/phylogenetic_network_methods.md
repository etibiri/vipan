# Phylogenetic network analysis — methods documentation

This document describes, in full, the phylogenetic-network method added to
the workflow (`scripts/build_phylogenetic_network.py`,
`scripts/classify_accessions.py`,
`scripts/analyze_network_vs_recombination.py`,
`workflow/phylogenetic_network.smk`). It is the methodological reference
cited from those scripts' docstrings and from the manuscript text in
`docs/manuscript_network_section_draft.md`. Every number quoted here is
reproduced from the workflow's own output files; none is invented.

## 1. Why a split network rather than a single tree

A single bifurcating tree (NJ, ML, or otherwise) forces every relationship
between genomes into a strictly hierarchical, non-reticulate structure. For
a dataset in which extensive recombination has already been documented by
independent methods (3SEQ: 6,196,942 recombinant triplets among 31,255,600
tested; RDP4: 126 robust events after accession expansion — see
`results/.../recombination/`), a tree is known a priori to be an
oversimplification for at least some genomes: a recombinant genome, by
definition, does not have a single coherent phylogenetic history, so
forcing it onto one branch of one tree discards exactly the signal the
recombination analyses were designed to detect.

A phylogenetic (split) network relaxes this constraint: instead of a single
set of nested bipartitions (a tree), it represents the data as a
possibly-larger set of bipartitions ("splits") of the taxon set, each with
a weight, drawn so that mutually incompatible splits appear as crossing
edges (reticulations) in the plotted network rather than being forced into
a single hierarchy. This is the standard rationale for using Neighbor-Net /
split-decomposition-family methods (Bandelt & Dress 1992, *Adv. Appl.
Math.* 13:92–105; Bryant & Moulton 2004, *Mol. Biol. Evol.* 21:255–265;
Dress & Huson 2004, *Appl. Bioinformatics* 3:7–8; Huson & Bryant 2006,
*Mol. Biol. Evol.* 23:254–267) on recombination-rich or otherwise
non-treelike datasets, and is why SplitsTree-family tools are routinely
used alongside, not instead of, recombination detection software in
begomovirus/geminivirus studies.

## 2. Why a from-scratch implementation, and what it is and is not

**Constraint.** This workflow's execution environments (the analysis
sandbox and the connected HPC-adjacent shell) have no outbound network
access to package registries (PyPI, CRAN, Ubuntu package mirrors all
return 403/blocked-by-allowlist errors when probed). Installing R +
`phangorn`/`ape`, or a SplitsTree binary, was therefore not possible in
this environment. Rather than skip the analysis, or claim to have run a
canonical NeighborNet implementation without being able to verify it, a
Neighbor-Net-family circular split network method was implemented from
scratch in `scripts/build_phylogenetic_network.py`, using only the
scientific Python stack already available (numpy, scipy, pandas,
matplotlib, networkx). This is stated explicitly here so the method's
scope and limitations are transparent to any reader or reviewer.

**What the implementation does** (see the script's own function-level
docstrings for exact algorithmic detail):

1. Compute a pairwise uncorrected p-distance matrix from the aligned
   FASTA (pairwise-deletion of gap/ambiguity columns — the same distance
   convention already used elsewhere in this workflow, e.g.
   `qc/distance_filter_decisions.tsv`).
2. Build a **circular taxon ordering** via classical neighbor-joining
   (Saitou & Nei 1987) followed by a depth-first traversal of the
   resulting tree. This substitutes for NeighborNet's own agglomerative
   ordering heuristic; it is a documented simplification, not the
   original NeighborNet ordering algorithm.
3. Enumerate **circular splits**: contiguous arcs of the circular
   ordering. The full circular split system for *n* taxa has
   *n*(*n*−1)/2 splits (≈54,000 for *n*≈329), which is computationally
   intractable for dense non-negative least squares (NNLS) fitting in
   this environment (a dense design matrix at that scale would require
   tens of gigabytes and, by empirical benchmarking during development,
   well over an hour per dataset). The split system is therefore bounded
   by `--max-split-side` (K=5 in this workflow's configuration,
   `config.yaml: phylogenetic_network.max_split_side`): only splits whose
   smaller side contains at most K taxa are candidates. This yields 1,645
   candidate splits for the 329-taxon dataset and 1,580 for the 316-taxon
   dataset — tractable (≈90–100 s per dataset) while still covering all
   splits that isolate small-to-moderate clades, which is where
   recombination-driven local reticulation is expected to manifest as
   conflicting signal against the deeper backbone.
4. Fit **non-negative split weights** to the candidate circular splits by
   a single NNLS solve against the full pairwise-distance vector (all
   *n*(*n*−1)/2 pairs), using `scipy.optimize.nnls`.
5. Compute an **equal-angle network layout**: taxon *i*'s position is the
   vector sum, over every split separating *i* from a fixed reference
   taxon, of (split weight) × (unit vector at the mean angular position of
   the side of the split *not* containing the reference taxon). This is
   the standard simplified equal-angle heuristic (cf. Dress & Huson 2004).
6. Plot the network: an outer skeleton connecting circularly-adjacent
   taxa, translucent chords for the largest-weight internal (non-trivial)
   splits, taxa coloured by their classification-table group (see §3).
7. Write NEXUS `SPLITS` and PHYLIP distance-matrix files
   (SplitsTree/Dendroscope-compatible), a splits table (`*.splits.tsv`,
   listing every retained split's weight and both sides' taxon lists),
   the circular order, and a JSON summary (`*.network_summary.json`)
   recording *n* taxa, alignment length, number of candidate/plotted/
   non-trivial splits, maximum split weight, and the NNLS sum-of-squared
   residual.

**What it is not.** This is *not* a validated reimplementation of the
canonical NeighborNet agglomeration algorithm (Bryant & Moulton 2004),
which was judged too complex to reimplement from scratch without a
reference implementation to cross-check against — attempting that risked
introducing silent, unvalidatable bugs. The substitution (NJ-derived
ordering in place of NeighborNet's own agglomeration step, followed by
the same downstream circular-split NNLS fitting and equal-angle layout
used by NeighborNet-family methods) is a documented approximation. It was
verified on a small (6-taxon) toy dataset with a known split structure
before being run at full scale (NNLS residual ≈5×10⁻³², i.e. an exact fit
on data that is exactly tree-like by construction, confirming the
NNLS/layout machinery is implemented correctly). At full scale, the
method achieves a fit statistic (§4) of 94–95%, comparable to typical
NeighborNet fit values reported in the literature for moderately
tree-like viral datasets, but this comparison should be read qualitatively
— not as a claim of numerical equivalence to SplitsTree/NeighborNet output
— because of the ordering-heuristic substitution and the bounded (K=5)
split system (point 3 above), the latter of which is a documented lower
bound on achievable fit, not an unbounded fit.

## 3. Classification of accessions (species / genetic group / lineage)

`scripts/classify_accessions.py` builds the accession → group table used
to colour and interpret both networks (full column-by-column derivation
is in that script's own module docstring). In summary, and using
**only** terminology the data actually support:

- **species** — the GenBank/NCBI submitted organism name, after two
  purely mechanical text normalizations (stripping a trailing
  bracket-style annotation; harmonizing the "Sweetpotato"/"Sweet potato"
  spelling variants used inconsistently across submitters for the same
  species). No species are merged, split, or reassigned by this step.
- **coarse_cluster** — average-linkage hierarchical clustering of the
  pairwise p-distance matrix, cut at 0.09 (≈91% identity), the value
  commonly cited as the species-demarcation threshold for begomoviruses
  (Brown et al. 2015, *Arch. Virol.* 160:1593–1619; Muhire et al. 2013,
  SDT). Computed independently per dataset (329 and 316).
- **species_or_group** — for each coarse cluster, the modal submitted
  species name, *provided* it covers ≥70% of the cluster's members;
  otherwise the cluster is reported as "Mixed cluster *N*" rather than
  being forced under a possibly-incorrect single species label. This is
  the deliberate no-invented-groups safeguard: a cluster is only ever
  named after a species when the data support it at a stated purity
  threshold.
- **lineage** — coarse clusters exceeding 25 members (only the dominant
  "Sweet potato leaf curl virus" cluster, in both datasets) are further
  partitioned by average-linkage clustering restricted to that cluster's
  own distance submatrix, into ⌈*n*/25⌉ locally-coherent sub-groups,
  labelled "*species* genetic group *k*". **The term "subspecies" is never
  used** — it is not a rank recognised for this genus by the ICTV. This
  step exists because the GenBank "Sweet potato leaf curl virus" label is
  known, from the coarse clustering itself, to span genetically
  heterogeneous material; "genetic group" is used deliberately as a
  purely descriptive, rank-free term for these internally-clustered
  sub-partitions.
- **Recombination status** — `rdp4_robust_recombinant` /
  `n_rdp4_robust_events_as_recombinant` and `strict_validation_tier` /
  `strict_validation_class` are copied unchanged from this workflow's own
  prior outputs (`recombination/rdp_events_robust.tsv`,
  `validation/event_validation_calls_strict.tsv`). No new recombination
  detection happens in this step; it only cross-references
  already-validated calls.

Resulting composition (from `accession_classification.tsv`, restricted to
`in_329_dataset` / `in_316_dataset` respectively):

| | 329-dataset | 316-dataset |
|---|---|---|
| n taxa | 329 | 316 |
| n species/groups (incl. "Mixed cluster *N*") | 16 | 15 |
| n lineages (species/groups + genetic-group sub-splits) | 22 | 21 |
| n RDP4-robust recombinant accessions present | 95 | 95 |

(The RDP4-robust recombinant count is identical between the two datasets
because none of the 13 genomes excluded by the distance/quality filter
between the 329- and 316-sequence datasets happened to carry an RDP4-robust
recombinant flag — a coincidence of which genomes failed QC, not a
methodological artifact.)

## 4. Fit statistic

A SplitsTree-style fit percentage is reported:

  fit% = 100 × (1 − SSR / Σ dᵢⱼ²)

where SSR is the NNLS sum-of-squared residual (from
`*.network_summary.json`) and Σ dᵢⱼ² is the sum of squared pairwise
p-distances (from `*.pairwise_distances.tsv`), both over the same
*n*(*n*−1)/2 pairs used in the NNLS fit. This is computed in
`scripts/analyze_network_vs_recombination.py::fit_statistic()`.

| | 329-dataset | 316-dataset |
|---|---|---|
| Alignment length (nt) | 3,823 | 3,757 |
| Candidate circular splits (K≤5) | 1,645 | 1,580 |
| Splits retained for plotting (≥0.2% max weight) | 574 | 555 |
| ...of which non-trivial (both sides >1 taxon) | 427 | 415 |
| **Fit statistic** | **95.18%** | **94.32%** |

A fit in the low-to-mid 90s indicates a dataset that is predominantly
tree-like — consistent with the bulk of sweepovirus diversification being
driven by vertical descent — with a real but bounded amount of conflicting
(reticulate-compatible) signal, which is exactly the qualitative picture
expected for a dataset in which the great majority of genomes are not
themselves recombinant (recall: 95/329 = 28.9% and 95/316 = 30.1% of
genomes carry an RDP4-robust recombination flag) but a substantial
minority are.

## 5. Cross-referencing the network against RDP4/3SEQ

`scripts/analyze_network_vs_recombination.py` (run via
`python scripts/analyze_network_vs_recombination.py <phylogenetic_network
output dir>`, writing `network_recombination_crossref.json`) asks, for
each dataset independently, whether the taxa placed on the "conflicting"
side of the network's largest splits are enriched for RDP4-robust
recombinants relative to the rest of the dataset. Method:

1. Take the 20 largest-weight **non-trivial** splits (both sides >1
   taxon) in each network — these are the splits with the strongest NNLS
   support, i.e. the most confidently inferred bipartitions.
2. For each split, take its smaller ("minority") side — the
   biologically interpretable candidate clade — and compute what fraction
   of its members share the same `lineage` label (from §3). A split whose
   minority side is ≥90% single-lineage is called **"clean"** (congruent
   with the independently-derived lineage partition, i.e. tree-like at
   that split); a split whose minority side mixes ≥2 lineages more evenly
   is called **"conflicting"** (a candidate reticulation).
3. Pool the union of taxa ever placed on a minority side across the 20
   largest splits, and test (two-sided Fisher's exact test, computed from
   scratch via the hypergeometric distribution — no scipy dependency, see
   `fisher_exact_greater()` in the script — to avoid adding an
   unavailable dependency) whether RDP4-robust recombinants are enriched
   on these minority sides versus the rest of the taxa.

Results (from `network_recombination_crossref.json`):

| | 329-dataset | 316-dataset |
|---|---|---|
| "Clean" splits among top 20 | 19/20 | 18/20 |
| "Conflicting" splits among top 20 | 1/20 | 2/20 |
| Taxa on any top-20 minority side | 39 | 69 |
| ...of which RDP4-robust recombinant | 17 (43.6%) | 34 (49.3%) |
| Baseline recombinant rate (whole dataset) | 28.9% | 30.1% |
| Fisher's exact odds ratio | 2.11 | 2.98 |
| Fisher's exact p-value (two-sided) | 0.038 | 1.7×10⁻⁴ |

**Interpretation.** The overwhelming majority (19/20 and 18/20) of the
strongest-supported splits in both networks are congruent with the
independently-derived lineage classification — the networks are, at their
core, tree-like, as expected for a dataset in which most genomes are not
recombinant. At the same time, RDP4-robust recombinants are
significantly over-represented among the taxa placed on the
minority/conflicting side of the strongest splits, in both datasets
independently. This is the expected signature if genuine recombination
is contributing measurable reticulate signal to the network on top of an
otherwise clonal backbone: it does **not** by itself prove that any
individual conflicting split corresponds to a specific detected
recombination event (the network and RDP4/3SEQ are independent methods
operating on different signals — whole-genome distance-based topology
versus site-pattern-based breakpoint detection), but the systematic,
statistically significant co-occurrence across both independently-built
networks is genuine, data-supported corroborating evidence.

The single genuinely "conflicting" split in the 329-dataset network
(weight 0.0351, the 6th-largest non-trivial split) placed four accessions
— KJ476507.1, MK052983.1, NC_021719.1, NC_038465.1 — split evenly between
two lineages ("Mixed cluster 4" and "Sweet potato leaf curl Henan virus").
Three of the four (KJ476507.1, NC_021719.1, NC_038465.1) are RDP4-robust
recombinants; KJ476507.1 and NC_038465.1 additionally carry an existing
strict-tier phylogenetic classification from this workflow's own
validation step (both Tier B; NC_038465.1 = class D, "expected switch,
insufficient branch support"; KJ476507.1 = class E, "parent switch
without expected roles"). The two conflicting splits in the 316-dataset
network show the same pattern (4/5 and 3/5 of minority-side taxa are
RDP4-robust recombinants respectively). These are reported as concrete,
traceable examples, not as new recombination detections — the events
themselves were already called by RDP4/3SEQ and validated (or found
unresolved) by the existing strict-tier pipeline; the network analysis
here shows that the same accessions independently produce anomalous,
tree-incongruent placement in a whole-genome distance-based network,
which is additional, methodologically independent support for treating
those specific accessions' evolutionary history as reticulate rather than
purely bifurcating.

## 6. Reproducibility

The full pipeline (`workflow/phylogenetic_network.smk`, three rules:
`distance_matrix` → `classify_accessions` → `phylogenetic_network`, plus
convenience target `phylogenetic_networks_all`) is Snakemake-integrated,
resource-annotated consistently with the existing SLURM profile
(`profiles/slurm/config.yaml`), and runnable either on that profile
(`snakemake --profile profiles/slurm phylogenetic_networks_all`) or
locally (`snakemake -c4 --use-conda phylogenetic_networks_all`). It must
be added to the main `Snakefile` via `include:
"workflow/phylogenetic_network.smk"` (matching this project's existing
opt-in pattern for `.smk` files — none are included by default).
Cross-referencing against RDP4/3SEQ (§5) is a separate, explicit script
invocation (`scripts/analyze_network_vs_recombination.py`) rather than a
Snakemake rule, since it is an interpretive/reporting step rather than a
pipeline stage producing a manuscript-ready artifact file; it can be
wired into the Snakemake DAG as an additional rule if continuous
re-generation on every workflow run is desired.
