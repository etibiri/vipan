# Workflow reproductible de recombinaison SPLCV

## Principe général

Le workflow ne modifie jamais l'analyse historique située dans `results/`.
Chaque reprise reçoit un `run_id` différent et écrit dans
`results/<run_id>/`. Le fichier `config.yaml` contient tous les seuils qui
influencent une décision biologique.

La chaîne complète est :

1. normalisation des identifiants FASTA et contrôle de longueur/ambiguïtés ;
2. orientation des génomes circulaires sur le motif `TAATATTAC` ;
3. déduplication exacte optionnelle ;
4. premier alignement MAFFT utilisé pour le diagnostic ;
5. calcul des distances nucléotidiques et exclusion auditée des séquences
   éloignées ou de mauvaise qualité ;
6. réalignement de novo des séquences retenues ;
7. résumé de qualité de l'alignement final ;
8. analyse autonome 3SEQ ;
9. import et filtrage du rapport RDP4 ;
10. concordance événementielle RDP/3SEQ ;
11. SimPlot/BootScan, phylogénies segmentaires et appels stricts ;
12. analyse des sites variables et parcimonie-informatifs ;
13. manifeste final contenant les empreintes SHA-256 des entrées et résultats.

## Pourquoi deux alignements ?

Le premier alignement sert uniquement à identifier les séquences qui ont une
couverture insuffisante, trop de gaps ou une distance excessive. Les séquences
retenues sont ensuite débarrassées des gaps et réalignées depuis les séquences
non alignées. Les coordonnées RDP et 3SEQ sont donc toujours définies sur
l'alignement final, jamais sur l'alignement diagnostique.

## Filtrage des séquences éloignées

Par défaut, l'ancre est le **médoïde** : la séquence ayant la plus faible
distance médiane envers le reste du jeu ayant passé le contrôle de qualité.
Une séquence est exclue si sa p-distance à ce médoïde dépasse
`distance_filter.max_p_distance` ou si le nombre de positions A/C/G/T
comparables est insuffisant.

Pour une analyse centrée sur une lignée prédéfinie, utiliser plutôt :

```yaml
distance_filter:
  mode: reference
  reference_accession: NC_XXXXXXXX.X
  max_p_distance: 0.20
```

Le seuil de 0,25 présent dans le fichier d'exemple n'est pas une constante
taxonomique universelle. Il doit être justifié par la distribution des
distances, la portée taxonomique de l'étude et, idéalement, une phylogénie
préliminaire. Un seuil trop strict peut supprimer des lignées parentales
réellement impliquées dans des recombinaisons interspécifiques.

Les sorties suivantes doivent être examinées avant de lancer RDP :

- `qc/distance_filter_decisions.tsv` ;
- `qc/pairwise_distances_initial.tsv` ;
- `qc/distance_filter.pdf` ;
- `qc/final_alignment_summary.tsv` ;
- `alignment/initial_alignment_excluded.fasta`.

Une liste d'accessions peut être imposée avec `manual_keep` ou
`manual_exclude`. Une inclusion manuelle peut neutraliser le seuil de distance,
mais ne neutralise pas les contrôles fondamentaux de longueur, gaps et
ambiguïtés.

## Préparation automatisée de l'alignement et analyse 3SEQ

Créer une fois l'environnement qui pilote Snakemake :

```bash
conda env create -f envs/workflow.yaml
conda activate vipan-recombination-workflow
```

Adapter d'abord `run_id`, le chemin du FASTA brut et les chemins de 3SEQ dans
`config.yaml`, puis effectuer un test à blanc :

```bash
snakemake --use-conda --cores 8 -n prepare_filtered_alignment
```

Produire d'abord l'alignement filtré et ses rapports, sans lancer 3SEQ :

```bash
snakemake --use-conda --cores 8 \
    --rerun-incomplete \
    --printshellcmds \
    prepare_filtered_alignment
```

Après inspection et validation du seuil, lancer 3SEQ :

```bash
snakemake --use-conda --cores 8 \
    --rerun-incomplete \
    --printshellcmds \
    prepare_recombination_inputs
```

Sur Slurm avec le plugin d'exécution installé :

```bash
snakemake --profile profiles/slurm prepare_recombination_inputs
```

Le fichier à importer dans RDP4 est :

```text
results/<run_id>/alignment/final_alignment.fasta
```

## Frontière RDP4

Dans la configuration par défaut, RDP4 GUI constitue le seul point manuel. Il
faut ouvrir **exactement** l'alignement final produit par le workflow, conserver
les paramètres RDP avec le manuscrit, puis exporter le rapport CSV vers le
chemin `rdp.report_csv` défini dans `config.yaml`.

Il ne faut jamais réutiliser le rapport RDP de l'analyse à 327 séquences après
avoir filtré ou réaligné les séquences : les accessions, les coordonnées de
cassure et les tests multiples ont changé.

Lorsqu'une exécution RDP4 en ligne de commande sera validée et disponible sur
le cluster, cette frontière pourra être remplacée par une règle dédiée. Cette
commande devra auparavant être validée sur un petit jeu témoin, car ses options
dépendent de la version installée.

## Analyse complète après export RDP

Une fois le rapport CSV placé au bon endroit :

```bash
snakemake --use-conda --cores 8 \
    --rerun-incomplete \
    --printshellcmds
```

Sur Slurm :

```bash
snakemake --profile profiles/slurm
```

Les résultats principaux seront écrits dans :

```text
results/<run_id>/recombination/
results/<run_id>/validation/
results/<run_id>/provenance/run_manifest.json
```

## Reprise avec un autre seuil

Pour comparer plusieurs règles d'exclusion, ne pas écraser une analyse
existante. Modifier simultanément le seuil et le `run_id` :

```yaml
run_id: splcv_refiltered_d020
distance_filter:
  mode: medoid
  max_p_distance: 0.20
```

On peut alors comparer les nombres de séquences, la longueur de l'alignement,
les sites informatifs et les appels de recombinaison entre les différents
dossiers de résultats.

## Contrôles recommandés avant publication

1. Vérifier visuellement toutes les séquences exclues et leur position dans une
   phylogénie préliminaire.
2. Réaliser au moins une analyse de sensibilité avec deux seuils de distance.
3. Déclarer le nombre de séquences exclues pour chaque motif.
4. Distinguer les exclusions de qualité des exclusions taxonomiques.
5. Archiver le FASTA brut, le FASTA final, `config.yaml`, les environnements
   Conda, le DAG Snakemake et le manifeste SHA-256.
6. Ne pas conclure qu'une séquence est « erronée » uniquement parce qu'elle est
   éloignée : elle peut représenter une lignée biologique légitime.

## Vérification du workflow

Compiler les scripts et lancer les tests unitaires :

```bash
python -m py_compile scripts/*.py
python -m unittest discover -s tests -p 'test_*.py'
```

Valider la syntaxe et le DAG Snakemake :

```bash
snakemake --lint
snakemake --use-conda --cores 1 -n prepare_recombination_inputs
```
rison.tsv
results/sensitivity_medoid_vs_NC004650/discordant_accessions.tsv
results/sensitivity_medoid_vs_NC004650/filtering_strategy_comparison.pdf
```

Les métriques principales sont la concordance brute des décisions, l'indice de
Jaccard des ensembles conservés et le kappa de Cohen. Les accessions
`medoid_only` ou `reference_only` doivent être inspectées taxonomiquement avant
de choisir l'analyse principale.

Pour une comparaison de sensibilité de toute la chaîne de recombinaison, RDP4
doit être relancé indépendamment sur chacun des deux alignements finaux. Il ne
faut pas utiliser le même export RDP dans les deux bras.

Après avoir terminé les deux workflows complets, comparer automatiquement les
classes strictes A--F :

```bash
snakemake \
    --snakefile workflow/validation_sensitivity.smk \
    --use-conda \
    --cores 1
```

Cette étape produit une table appariée par accession, le nombre de changements
de classe et l'indice de Jaccard des ensembles soutenus par les classes A--C.

## Contrôles recommandés avant publication

1. Vérifier visuellement toutes les séquences exclues et leur position dans une
   phylogénie préliminaire.
2. Réaliser au moins une analyse de sensibilité avec deux seuils de distance.
3. Déclarer le nombre de séquences exclues pour chaque motif.
4. Distinguer les exclusions de qualité des exclusions taxonomiques.
5. Archiver le FASTA brut, le FASTA final, `config.yaml`, les environnements
   Conda, le DAG Snakemake et le manifeste SHA-256.
6. Ne pas conclure qu'une séquence est « erronée » uniquement parce qu'elle est
   éloignée : elle peut représenter une lignée biologique légitime.

## Vérification du workflow

Compiler les scripts et lancer les tests unitaires :

```bash
python -m py_compile scripts/*.py
python -m unittest discover -s tests -p 'test_*.py'
```

Valider la syntaxe et le DAG Snakemake :

```bash
snakemake --lint
snakemake --use-conda --cores 1 -n prepare_recombination_inputs
```
