# Validation ciblée des événements RDP/3SEQ

## Objectif

Le script `scripts/validate_recombination_events.py` transforme les événements
prioritaires du tableau de concordance en trois niveaux de preuve complémentaires :

1. **SimPlot-style** : identité nucléotidique glissante du recombinant envers
   chaque parent candidat ;
2. **BootScan-style** : fréquence à laquelle chaque parent est le terminal le
   plus proche du recombinant dans des arbres NJ construits après rééchantillonnage
   des colonnes de chaque fenêtre ;
3. **phylogénies segmentaires** : arbres de maximum de vraisemblance IQ-TREE 2
   pour le segment entre les cassures et son complément circulaire.

Le script déduit automatiquement `Tier_A` à `Tier_E` si la colonne
`Evidence_tier` n'existe pas. Avec les résultats actuels, `Tier_B` correspond
aux événements à concordance de cassures modérée/partielle ayant au moins un
parent concordant : c'est le premier lot à vérifier.

## Installation

```bash
conda env create -f envs/recombination_validation.yaml
conda activate vipan-recombination-validation
```

## Test court sur un événement

```bash
python scripts/validate_recombination_events.py \
    --events results/recombination/concordance_3seqraw_rdp_327.event_pairs.tsv \
    --alignment results/alignment/splcv_aligned.fasta \
    --outdir results/recombination/validation_test \
    --tiers Tier_B \
    --topology circular \
    --window 200 \
    --step 20 \
    --bootstrap-replicates 50 \
    --processes 4 \
    --context-size 8 \
    --max-events 1
```

Cette commande vérifie la sélection, les profils et les FASTA segmentaires
sans lancer IQ-TREE. Inspecter le PDF `simplot_bootscan.pdf` produit dans le
dossier de l'événement.

## Analyse de validation complète

```bash
srun --cpus-per-task=8 --mem=16G --time=12:00:00 --pty bash

conda activate vipan-recombination-validation

python scripts/validate_recombination_events.py \
    --events results/recombination/concordance_3seqraw_rdp_327.event_pairs.tsv \
    --alignment results/alignment/splcv_aligned.fasta \
    --outdir results/recombination/validation \
    --tiers Tier_B \
    --topology circular \
    --window 200 \
    --step 20 \
    --bootstrap-replicates 500 \
    --processes 8 \
    --rdp-candidates-per-role 2 \
    --context-size 8 \
    --run-iqtree \
    --iqtree-model MFP \
    --ufboot 1000 \
    --alrt 1000 \
    --iqtree-threads 8 \
    > results/logs/validate_recombination_events.log 2>&1
```

Sur un cluster géré par Snakemake, ajouter au `Snakefile` :

```python
include: "workflow/recombination_validation.smk"
```

Puis ajouter éventuellement dans `config.yaml` :

```yaml
recombination_validation:
  events: results/recombination/concordance_3seqraw_rdp_327.event_pairs.tsv
  alignment: results/alignment/splcv_aligned.fasta
  tiers: Tier_B
  topology: circular
  window: 200
  step: 20
  bootstrap_replicates: 500
  rdp_candidates_per_role: 2
  context_size: 8
  run_iqtree: true
  reuse_existing: false
  formats: png,pdf
  threads: 8
  mem_mb: 16000
  runtime_minutes: 720
```

Lancement :

```bash
snakemake --use-conda --cores 8 validate_recombination_events -p
```

## Sorties

À la racine du dossier de sortie :

- `validation_manifest.tsv` : événements traités, parents et séquences de
  contexte sélectionnés ;
- `segment_tree_summary.tsv` : arbre, parent sélectionné le plus proche,
  distance patristique et support du clade recombinant-parent ;
- `event_validation_calls.tsv` : appel synthétique du changement de parent
  entre les deux segments ;
- un sous-dossier par événement.

Dans chaque sous-dossier :

- `simplot_input.fasta` : recombinant, parents et contextes ; ce fichier peut
  aussi être ouvert dans SimPlot++ pour une vérification graphique indépendante ;
- `selected_parents.tsv` ;
- `simplot_data.tsv` et `bootscan_data.tsv` ;
- `simplot_bootscan.pdf` et `simplot_bootscan.png` ;
- `segment_trees/RDP/` et `segment_trees/3SEQ/`, contenant les alignements
  segmentaires, les fichiers IQ-TREE et les figures d'arbre.

## Lecture scientifique

Pour le schéma de cassures RDP, le segment A (entre les cassures) doit se
rapprocher du parent mineur et le segment B du parent majeur. Pour le schéma
3SEQ, les parents 1 et 2 peuvent permuter entre A et B, car leur orientation
mineur/majeur n'est pas imposée.

`Supported_expected_parent_switch` exige :

- deux parents les plus proches différents entre les segments ;
- les rôles parentaux attendus ;
- SH-aLRT au moins 80 et UFBoot au moins 95 sur les deux clades concernés.

Les arbres IQ-TREE sont non enracinés. Les supports recombinant-parent sont
donc extraits de la plus petite bipartition soutenue qui regroupe les deux
taxons, en examinant les deux côtés de chaque branche. Cette procédure évite de
dépendre de la racine arbitraire utilisée pour écrire le fichier Newick.

Cet appel est une aide au tri, pas une preuve automatique. Un événement sera
retenu comme bien validé seulement si le changement observé est compatible
entre le profil SimPlot, le BootScan, les arbres segmentaires et les cassures,
sans séquence de contexte plus proche expliquant le signal. Les arbres de
segments très courts, peu variables ou fortement lacunaires doivent être
interprétés comme non concluants.

Par défaut, les séquences de contexte sont les huit génomes complets les plus
proches du recombinant, hors parents retenus. Pour imposer un échantillonnage
taxonomique contrôlé, créer un fichier contenant une accession par ligne et
ajouter `--context-accessions chemin/contextes.txt`. Les accessions absentes de
l'alignement commun sont ignorées.

Après une mise à jour du script, les arbres et profils déjà produits peuvent
être réanalysés sans recalculer IQ-TREE ni le BootScan en ajoutant
`--reuse-existing` et en conservant exactement le même dossier `--outdir`.
