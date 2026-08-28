# Draft manuscript text — phylogenetic network analysis

**How to use this file.** These paragraphs are written to slot into the
existing manuscript structure (`manuscript_v3.md` / the corresponding
`.docx`), which is numbered 2.1–2.11 (Materials and methods), 3.1–3.6
(Results), 4.1–4.7 (Discussion). Suggested insertion points are given
before each block. All figures/tables referenced (`Figure 6`, `Figure 7`,
`Table 5`) are placeholders — renumber to fit the manuscript's final
figure/table sequence. Every number below is reproduced from
`results/splcv_refiltered_medoid_deterministic_v2/phylogenetic_network/`
(`*.network_summary.json`, `network_recombination_crossref.json`,
`accession_classification.tsv`); none is invented. Full methodological
detail and citations are in `docs/phylogenetic_network_methods.md`. Note
that this manuscript already correctly refers to "RDP4" throughout
(Section 2.5, "RDP4 screening and robust-event filtering") — the Abstract
and Keywords in the current draft still read "RDP5" and should be
corrected to "RDP4" when the live document is next editable (see the
project-wide RDP5→RDP4 correction note delivered separately).

---

## Materials and methods — insert as new subsection 2.12 (after 2.11 Workflow reproducibility)

### 2.12. Phylogenetic network analysis

Because a fraction of the dataset is recombinant by construction (95 of
329 genomes, 28.9%, carried at least one RDP4-robust recombination event;
§2.5), a strictly bifurcating tree was judged an inadequate summary of
genome-wide relationships for this dataset: forcing a recombinant genome
onto a single tree branch is expected to distort or discard exactly the
mosaic signal that 3SEQ and RDP4 were used to detect. We therefore
constructed circular split networks (the Neighbor-Net family;
[Bryant & Moulton, 2004]; [Bandelt & Dress, 1992]), which represent
genome-wide relationships as a set of weighted, potentially incompatible
splits rather than a single hierarchy, so that unresolved or reticulate
signal appears as crossing edges rather than being forced into an
arbitrary bifurcation.

Networks were built independently for the 329-genome dataset (before
distance-based filtering) and the 316-genome dataset (the filtered set
used for the primary RDP4/3SEQ analysis; §2.3), so that the effect of
the 13 excluded genomes on overall network structure could also be
assessed. Because the standard network-construction toolchain
(SplitsTree; R/phangorn) could not be installed in the execution
environment used for this analysis (no outbound access to package
repositories), the method was implemented directly in Python using the
existing scientific stack (numpy/scipy/pandas), following the standard
circular-split-network construction pipeline: (i) a pairwise
uncorrected-p-distance matrix (pairwise deletion of gaps/ambiguous
sites, consistent with the distance convention used in §2.3); (ii) a
circular taxon ordering derived from a neighbor-joining tree
([Saitou & Nei, 1987]) — used here as a documented substitute for
NeighborNet's own internal ordering heuristic; (iii) enumeration of
circular splits (contiguous arcs of the ordering), bounded to splits
whose smaller side contains at most five taxa for computational
tractability (the unbounded circular split system has *n*(*n*−1)/2
candidate splits, intractable for dense non-negative least-squares
fitting at *n*≈329 in this environment); (iv) non-negative
least-squares fitting of split weights to the full pairwise-distance
matrix (`scipy.optimize.nnls`); and (v) an equal-angle network layout.
Outputs — NEXUS split blocks, PHYLIP distance matrices, split tables,
and publication-resolution PDF/SVG/PNG figures — were generated with
`scripts/build_phylogenetic_network.py` and are reproducible end-to-end
via the Snakemake rules in `workflow/phylogenetic_network.smk`,
including on the project's existing SLURM cluster profile. The method,
its parameters, and its explicitly documented limitations relative to a
canonical NeighborNet implementation are described in full in
`docs/phylogenetic_network_methods.md`.

Each genome was independently assigned a species/genetic-group label
using only data-driven or unmodified-metadata evidence, to avoid
asserting taxonomic distinctions the data do not support: (i) the
GenBank-submitted organism name, after purely mechanical text
normalization; (ii) an average-linkage hierarchical clustering of the
pairwise p-distance matrix cut at 0.09 substitutions/site (≈91%
identity, the demarcation value commonly used for begomovirus species;
[Brown et al., 2015]; [Muhire et al., 2013]), computed independently for
each dataset; a cluster was labelled with its modal submitted species
name only if that name covered ≥70% of the cluster's members, otherwise
it was reported as an unlabelled "mixed cluster" rather than assigned a
potentially incorrect species; and (iii), for the one oversized cluster
in both datasets (corresponding to the GenBank label "Sweet potato leaf
curl virus", already known to be genetically heterogeneous), a further
partition into internally coherent "genetic groups" by clustering
restricted to that cluster's own distance submatrix. We deliberately use
"genetic group" and "lineage" rather than "subspecies" throughout, as
subspecies is not a taxonomic rank recognized for this genus by the
ICTV. This classification (`scripts/classify_accessions.py`) also
carries forward, unchanged, each accession's existing RDP4-robust
recombinant status and strict-tier phylogenetic validation class from
§§2.5–2.9, allowing direct cross-referencing between network topology
and the recombination-detection results without any new recombination
calling.

To test whether network structure was compatible with the independently
detected recombination signal, for each dataset we took the 20
largest-weight non-trivial splits (i.e., the most strongly supported
candidate bipartitions) and classified each as "clean" if its
minority-side taxa shared a single genetic-group/lineage label at ≥90%
purity, or "conflicting" otherwise. We then tested, by two-sided Fisher's
exact test, whether RDP4-robust recombinant accessions were
over-represented among the pooled minority-side taxa of these 20 splits
relative to the rest of the dataset (`scripts/analyze_network_vs_
recombination.py`).

---

## Results — insert as new subsection 3.7 (after 3.6 Alignment determinism...)

### 3.7. Circular split networks are predominantly tree-like, with recombination-enriched local reticulation

The 329- and 316-genome circular split networks (Figure 6A–B) each
retained several hundred splits above the reporting threshold (574 and
555 splits respectively, of which 427 and 415 were non-trivial,
i.e. separated more than one taxon from the rest) from a bounded
candidate set of 1,645 and 1,580 circular splits. Both networks achieved
a high overall fit to the underlying pairwise-distance data (95.2% for
the 329-genome network; 94.3% for the 316-genome network; fit defined as
100 × [1 − residual sum of squares / total sum of squared distances]),
indicating that genome-wide relationships in this dataset are
predominantly, but not purely, tree-like.

Independently derived species/genetic-group clustering (§2.12) resolved
16 species-or-group labels and 22 genetic-group-level lineages in the
329-genome dataset, and 15 and 21 respectively in the 316-genome dataset
(Table 5; full accession-level assignments in Supplementary Table
S[accession_classification]). Among the 20 largest-weight (most strongly
supported) splits in each network, the great majority were congruent
with this independently derived classification: 19/20 in the
329-genome network and 18/20 in the 316-genome network placed a
lineage-pure clade (≥90% single-lineage purity) on the split's minority
side, consistent with the high overall fit statistic and with a
predominantly clonal/vertical backbone.

The remaining splits — one in the 329-genome network, two in the
316-genome network — mixed two or more genetic groups on their minority
side and are the clearest network-level candidates for reticulate,
recombination-compatible signal. Pooling the minority-side taxa across
all 20 top splits in each network and cross-referencing against the
RDP4-robust recombination calls from §2.5 showed a statistically
significant enrichment for recombinant accessions in both datasets
independently: 17/39 minority-side taxa (43.6%) were RDP4-robust
recombinants in the 329-genome network, versus a 28.9% baseline rate
across the whole dataset (Fisher's exact odds ratio = 2.11, two-sided
*p* = 0.038); 34/69 (49.3%) in the 316-genome network, versus a 30.1%
baseline (odds ratio = 2.98, *p* = 1.7 × 10⁻⁴).

The single conflicting split in the 329-genome network (relative weight
0.035, the sixth-largest non-trivial split) placed four accessions —
KJ476507.1, MK052983.1, NC_021719.1, and NC_038465.1 — split evenly
between two genetic groups ("Mixed cluster 4" and "Sweet potato leaf
curl Henan virus"). Three of these four (KJ476507.1, NC_021719.1,
NC_038465.1) are RDP4-robust recombinants; KJ476507.1 and NC_038465.1
also carry an existing strict-tier phylogenetic classification from
§2.9 (both Tier B; NC_038465.1 in class D, "expected switch, insufficient
branch support"; KJ476507.1 in class E, "parent switch without expected
roles"). The two conflicting splits in the 316-genome network show the
same pattern, each placing 3–4 of 5 minority-side taxa among the
RDP4-robust recombinants.

---

## Discussion — insert as new subsection 4.4bis, or fold into 4.4 "Relationship to known sweepovirus evolution"

Suggested placement: immediately after the existing §4.4 paragraph(s), as
a new paragraph (or, if the editors prefer a single combined
subsection, retitle 4.4 to "Relationship to known sweepovirus evolution
and network structure").

The phylogenetic network analysis provides an independent, whole-genome,
distance-based perspective on the same dataset examined by 3SEQ, RDP4,
and the segment-level phylogenetic validation reported above (§§2.4–2.10,
3.2–3.5). The two are methodologically orthogonal: RDP4 and 3SEQ detect
localized breakpoint/mosaicism signal from site patterns, and the strict
phylogenetic validation tests specific candidate breakpoints against
segment-level maximum-likelihood trees, whereas the split network
summarizes genome-wide pairwise distance structure without reference to
any specific breakpoint hypothesis. That RDP4-robust recombinants are
nonetheless significantly over-represented among the taxa producing
network-level topological conflict, in both the 329- and 316-genome
datasets independently (odds ratios 2.11–2.98, both *p* < 0.05), is
therefore corroborating rather than circular evidence: it shows that
genomes flagged as recombinant by site-pattern-based methods also tend
to occupy anomalous, tree-incongruent positions in a network built from
an entirely different signal (overall pairwise distance).

At the same time, the networks were overwhelmingly tree-like at their
best-supported splits (19/20 and 18/20 of the twenty largest-weight
splits were congruent with independently derived genetic-group
boundaries; overall fit 94–95%), and only one and two splits,
respectively, showed clear topological conflict. This is consistent with
the picture built from the concordance and strict-validation analyses
above (§§3.3–3.4): most of the genome-wide phylogenetic signal in this
dataset reflects vertical descent within well-resolved genetic groups,
and only a minority of genomes — a similar order of magnitude to the 21
priority recombinants examined by strict phylogenetic classification, of
which two received full support (§3.4) — contribute detectable
reticulate signal at the whole-genome network level. The network
analysis does not, by itself, confirm or refute any individual
breakpoint call; rather, it offers a second, methodologically
independent line of evidence that the small set of genomes already
flagged by RDP4/3SEQ as recombinant are also, as a group, associated with
measurable, statistically significant network-level conflict, which is
the pattern expected if the detected recombination signal reflects real
underlying reticulate evolutionary history rather than being an artifact
specific to one detection method.

We note that this network-based cross-check, by design, cannot and does
not identify recombination events that neither 3SEQ nor RDP4 already
flagged: the enrichment test is defined relative to the existing
RDP4-robust call set (§2.5), not as an independent recombination-calling
method. Its contribution is corroborative and integrative — connecting a
whole-genome structural signature to breakpoint-level detections already
made — rather than a fourth, independent detection method in its own
right.

---

## Suggested new items for the manuscript's methods/results figure and table list

- **Figure 6.** Circular split networks for the 329-genome (panel A) and
  316-genome (panel B) datasets, taxa coloured by genetic group/lineage
  (source: `.../phylogenetic_network/{329,316}_dataset/splcv_{329,316}.network.pdf`).
- **Table 5.** Accession-level species/genetic-group/lineage
  classification and cross-reference to RDP4-robust recombinant status
  and strict-validation tier (source: `accession_classification.tsv`,
  reduce/rename columns as needed for the manuscript's supplementary
  table conventions).
- **Optional supplementary table.** Full network fit statistics and
  clean/conflicting split counts per dataset (source:
  `network_recombination_crossref.json`), if reviewers request
  quantitative support for the network-vs-RDP4 cross-referencing beyond
  what is stated in text.
