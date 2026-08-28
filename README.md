# vipan — Sweepovirus recombination validation workflow

Reproducible, auditable Snakemake workflow for detecting, filtering, and
phylogenetically validating recombination events in sweet
potato-infecting begomoviruses (sweepoviruses), and for building
circular split (Neighbor-Net-family) phylogenetic networks of the same
genome collection. This is the code and data repository accompanying the
manuscript *"Event-level concordance and phylogenetic validation reveal a
wide gap between recombination detection and confirmation in sweet
potato-infecting begomoviruses"*.

**Repository:** https://github.com/etibiri/vipan

---

## 1. What this workflow does

1. Standardizes and quality-filters a public collection of complete
   sweepovirus genomes (circular-origin reorientation, deterministic
   MAFFT alignment, pairwise-distance-based filtering).
2. Screens for recombination independently with standalone **3SEQ** and
   **RDP4** (RDP4 is a GUI application and is the workflow's one
   deliberate manual boundary — see §5).
3. Matches 3SEQ and RDP4 events one-to-one by breakpoint and
   recombinant-accession concordance, and classifies match quality.
4. Phylogenetically validates priority (Tier A/B) candidate events using
   sliding-window similarity, bootstrap scanning, and segment-level
   maximum-likelihood trees for the recombinant region and its circular
   complement.
5. Builds circular split (Neighbor-Net-family) phylogenetic networks of
   the full genome collection, classifies every accession into
   species/genetic-group/lineage labels, and cross-references network
   topology against the RDP4/3SEQ recombination calls.
6. Produces every figure and table used in the manuscript from raw
   inputs, with a full, checksum-verifiable provenance record.

## 2. Prerequisites

- Linux or macOS (developed and tested on Linux; the SLURM profile
  targets an HPC cluster).
- [conda](https://docs.conda.io) or [mamba](https://mamba.readthedocs.io)
  (mamba strongly recommended — it is what `envs/workflow.yaml`
  installs and what the environments below assume).
- [Snakemake](https://snakemake.readthedocs.io) ≥ 8 (installed via
  `envs/workflow.yaml`, below).
- **RDP4** ([Martin et al. 2021](https://doi.org/10.1093/ve/veaa087)) —
  Windows GUI application (runs under Wine on Linux/macOS). This is the
  workflow's only non-command-line, manual step; see §5.
- **3SEQ** binary and p-value table (`resources/3seq/`, not included in
  this repository — see [the 3SEQ website](https://mol.ax/software/3seq/)
  for the academic-use license and download; place the binary at
  `resources/3seq/3seq` and the table at
  `resources/3seq/pvalue_tables/pvaluetable_700`, matching `config.yaml:
  threeseq`).
- No internet access is required to *run* the workflow once genomes and
  3SEQ/RDP4 are in place; the accession list itself is derived from
  public GenBank records (`data/raw/splcv_all-0726.tsv` /
  `.fasta`, already included in this repository).

## 3. Project structure

```
Snakefile                    Main Snakemake entry point (default target: `all`)
config.yaml                  Default run configuration
configs/                     Config overlays (--configfile), e.g. medoid vs.
                              fixed-reference sensitivity arms
workflow/*.smk               Optional rule modules (recombination_validation,
                              filtering_sensitivity, validation_sensitivity,
                              phylogenetic_network) -- NOT included by default;
                              see each file's docstring for the one-line
                              `include:` needed to wire it into the DAG
envs/*.yaml                  Conda environments (one per `conda:` directive)
profiles/slurm/config.yaml   SLURM executor profile for `--profile profiles/slurm`
scripts/                     All Python/shell scripts invoked by the rules
data/raw/                    Input genome FASTA + accession metadata TSV
data/rdp/                    RDP4 GUI CSV exports (the one manual-boundary input)
results/<run_id>/            All pipeline outputs, one subtree per run_id
docs/                        Methodology write-ups
figures_and_code/            Manuscript figure-generation script and outputs
zenodo_release/               Curated, self-contained subset of this repository
                              prepared for Zenodo archival (see its own README)
tests/                       Workflow reproducibility tests
```

Every biologically distinct run uses a new `run_id` (set in `config.yaml`
or a `--configfile` overlay); outputs are namespaced under
`results/<run_id>/` so re-running with different settings never
overwrites a prior run.

The manuscript's primary analysis is `run_id: splcv_refiltered_medoid_
deterministic_v2` (329 genomes in, 316 retained after distance
filtering), produced via the `configs/medoid_sensitivity.yaml` overlay
described below.

## 4. Environment installation

```bash
git clone https://github.com/etibiri/vipan.git
cd vipan

# Workflow driver (Snakemake + SLURM executor plugin)
mamba env create -f envs/workflow.yaml
mamba activate vipan-recombination-workflow

# The two analysis environments (recombination_core, recombination_validation)
# are created automatically per-rule by Snakemake's --use-conda; you do not
# need to create them manually unless running a script standalone outside
# Snakemake.
```

## 5. Data preparation

`data/raw/splcv_all-0726.fasta` and `.tsv` (complete-genome sequences and
their GenBank-derived metadata) are already included in this repository.
To rebuild them from scratch, or to re-derive an updated accession list,
see `scripts/fetch_sweepo.py`.

**RDP4 (the manual boundary).** Run the target `prepare_recombination_
inputs` first; it produces the final alignment and the standalone 3SEQ
result without requiring RDP4:

```bash
snakemake --configfile config.yaml configs/medoid_sensitivity.yaml \
    -c8 --use-conda prepare_recombination_inputs
```

Import `results/<run_id>/alignment/final_alignment.fasta` into RDP4,
run the standard method battery, and export the CSV report to the path
configured under `rdp.report_csv` in `config.yaml`
(`data/rdp/<run_id>_report.csv`). Then run the default target as below.

## 6. Running the workflow

**Primary manuscript run** (329→316 genomes, medoid-based distance
filtering):

```bash
snakemake --configfile config.yaml configs/medoid_sensitivity.yaml \
    -c8 --use-conda all
```

**On the project's SLURM cluster profile:**

```bash
snakemake --configfile config.yaml configs/medoid_sensitivity.yaml \
    --profile profiles/slurm --use-conda all
```

**Reference-sensitivity arm** (fixed reference NC_004650.1 instead of
the data-derived medoid — reproduces the manuscript's sensitivity
analysis, §2.3):

```bash
snakemake --configfile config.yaml configs/reference_NC_004650_sensitivity.yaml \
    -c8 --use-conda all
```

**Recombination validation, filtering-sensitivity, and phylogenetic
network analyses** are each defined in their own `workflow/*.smk` module
and are opt-in (not part of the default `all` target). Add the relevant
line to the `Snakefile` (near the top, after `configfile:`):

```python
include: "workflow/recombination_validation.smk"     # Tier A/B phylogenetic validation
include: "workflow/filtering_sensitivity.smk"         # medoid-vs-reference sensitivity comparison
include: "workflow/validation_sensitivity.smk"        # validation-strategy sensitivity comparison
include: "workflow/phylogenetic_network.smk"          # circular split networks (329- and 316-genome)
```

then run the corresponding target, e.g.:

```bash
snakemake --configfile config.yaml configs/medoid_sensitivity.yaml \
    -c8 --use-conda phylogenetic_networks_all
```

(or `--profile profiles/slurm` on the cluster). Each module's docstring
documents its own targets and inputs/outputs in full.

## 7. Phylogenetic network analysis

Circular split (Neighbor-Net-family) networks are built independently
for the 329-genome (pre-distance-filter) and 316-genome (filtered)
datasets by `scripts/build_phylogenetic_network.py`, orchestrated by
`workflow/phylogenetic_network.smk` (3 rules: `distance_matrix` →
`classify_accessions` → `phylogenetic_network`, plus the convenience
target `phylogenetic_networks_all`). Because standard network software
(SplitsTree, R/phangorn) could not be installed in the original analysis
environment, the method is implemented from scratch in Python
(numpy/scipy) with every simplification relative to the canonical
NeighborNet algorithm explicitly documented; see
`docs/phylogenetic_network_methods.md` for the full method, its
justification, and its limitations, and
`docs/manuscript_network_section_draft.md` for manuscript-ready
Materials & Methods / Results / Discussion text built from this
workflow's own output.

Cross-referencing the resulting networks against the RDP4/3SEQ
recombination calls is a separate script:

```bash
python scripts/analyze_network_vs_recombination.py \
    results/<run_id>/phylogenetic_network
```

which writes `network_recombination_crossref.json` alongside the network
outputs.

## 8. Figures and tables

All manuscript figures and summary tables are regenerated from workflow
outputs by:

```bash
cd figures_and_code
python generate_figures_and_tables.py
```

Edit the `BASE`/`OUT` paths at the top of that script to point at your
own `results/<run_id>/` tree if it differs from the default.

## 9. Where to find the main results

| Result | Path |
|---|---|
| Final filtered alignment (316 genomes) | `results/<run_id>/alignment/final_alignment.fasta` |
| Standalone 3SEQ recombination output | `results/<run_id>/recombination/3seq_output.rec` |
| RDP4-robust events (filtered) | `results/<run_id>/recombination/rdp_events_robust.tsv` |
| 3SEQ/RDP4 event concordance | `results/<run_id>/recombination/concordance.summary.tsv` |
| Strict phylogenetic validation calls (Tier A/B) | `results/<run_id>/validation/event_validation_calls_strict.tsv` |
| Per-event SimPlot/BootScan/segment-tree evidence | `results/<run_id>/validation/Tier_*__*/` |
| Accession → species/genetic-group/lineage classification | `results/<run_id>/phylogenetic_network/accession_classification.tsv` |
| Phylogenetic networks (329- and 316-genome) | `results/<run_id>/phylogenetic_network/{329,316}_dataset/*.network.{pdf,svg,png}` |
| Network-vs-recombination cross-reference | `results/<run_id>/phylogenetic_network/network_recombination_crossref.json` |
| Full run provenance (checksums, software versions) | `results/<run_id>/provenance/run_manifest.json` |
| Manuscript figures/tables | `figures_and_code/figure*.png`, `figures_and_code/table*.csv` |

A curated, self-contained subset of the above, prepared for Zenodo
archival, is in `zenodo_release/` (see its own `README.md`).

## 10. Citation

If you use this workflow, please cite the manuscript (full citation to
be added on publication) and/or the software release itself via
`CITATION.cff` / the archived Zenodo record (DOI to be minted on first
tagged GitHub release — see `zenodo_release/.zenodo.json`).

## 11. License

See `LICENSE`.

## 12. A note on nomenclature

RDP**4** (not RDP5) was the software actually used for recombination
detection in this workflow (`rdp.report_csv` inputs, §2.5 of the
manuscript). 
