# Jeu de données manioc : téléchargement et détection

Une seule commande sur la machine Ubuntu (accès Internet requis pour NCBI) :

```bash
cd ~/vipan_transfer
sudo apt install -y mafft        # si absent
export NCBI_EMAIL=votre.adresse@exemple.org   # demandé par NCBI, jamais écrit dans les fichiers
./cassava/run_cassava.sh A       # ADN-A ; relancer avec B pour l'ADN-B (analyses séparées)
```

Sorties dans `work/cassava/` :

- `raw/manifest.tsv` : accession, espèce, isolat, pays, longueur, segment (A, B ou unassigned), titre. À me renvoyer en premier : il sert à établir la liste des recombinants connus avec leurs sources.
- `segA/prepare_log.tsv` : pourquoi chaque séquence est gardée, écartée ou mise de côté (mauvais segment, longueur, doublon exact, motif d'origine introuvable, brin incertain).
- `segA/thinning_log.tsv` : séquences retirées après alignement (trop de trous, quasi identiques, plafond de séquences).
- `segA/alignment_final.fasta` : alignement utilisé pour la détection.
- `segA/detect/events.tsv` : événements 3SEQ et GENECONV au même format ; `timings.tsv` : temps de calcul.

Variables utiles : `MAX_SEQS` (120 par défaut, 0 = pas de plafond ; GENECONV devient lent au-delà), `THREADS`, `MIN_LEN` et `MAX_LEN` (2500 et 2900 par défaut).

## Ce qui est vérifié et ce qui ne l'est pas

Vérifié sur des données simulées : orientation des génomes circulaires (y compris quand la boucle d'origine est palindromique et que le motif apparaît sur les deux brins), rotation, retrait des doublons, alignement MAFFT, détection, format de sortie, et lecture d'un résumé NCBI factice.

Non vérifié : le téléchargement lui-même, faute d'accès réseau chez moi. Si une requête NCBI échoue, le message d'erreur s'affiche et le script réessaie six fois avant d'abandonner.

## Limites à connaître

- Le segment (A ou B) est déduit du titre de la notice GenBank. Les notices sans indication sont classées `unassigned` et exclues par défaut.
- Le brin est choisi par k-mers partagés, pas par le motif : à vérifier si `prepare_log.tsv` signale beaucoup de `strand_uncertain`.
- 3SEQ désigne un recombinant et deux parents ; GENECONV signale seulement un fragment partagé entre deux séquences, sans désigner le receveur.
- Les points de rupture 3SEQ sont les milieux des intervalles d'incertitude, en colonnes d'alignement.
- Les « recombinants connus » ne sont pas encore codés : ils viennent de la littérature et doivent être vérifiés accession par accession avant de servir de vérité terrain.
