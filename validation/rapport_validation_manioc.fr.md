# Validation de vipan sur les begomovirus de la mosaïque du manioc

Ce document rend compte de l'étape 6 du plan de refonte : mesurer ce que valent les détecteurs de recombinaison, et la couche de concordance de vipan, sur un jeu réel dont les événements ont déjà été publiés. Les simulations du banc v0 restent le bras où la vérité est exacte. Le manioc apporte autre chose : des génomes réels, une diversité qu'aucun simulateur ne reproduit, et des recombinants que la communauté reconnaît.

## Le jeu de données

1 842 notices GenBank de begomovirus infectant le manioc ont été récupérées, dont 1 032 ADN-A, 658 ADN-B et 152 sans segment indiqué dans le titre. Après contrôle de longueur, recadrage sur le nonanucléotide d'origine et retrait des doublons exacts, 974 ADN-A et 601 ADN-B ont été alignés avec MAFFT. Les deux segments sont traités séparément : un begomovirus bipartite n'a pas un génome, il en a deux, et les confondre dans un même alignement n'aurait aucun sens.

Le filtre post-alignement a dû être corrigé avant d'aller plus loin. Il comparait la proportion de tirets de chaque séquence à la largeur totale de l'alignement, or MAFFT avait ouvert 3 808 colonnes pour des génomes de 2 800 bases — près de 920 colonnes d'insertions portées par une poignée de séquences divergentes. Chaque génome affichait donc 26 % de tirets et tombait sous le seuil de 20 %. L'alignement entier était jeté. En retirant d'abord les colonnes occupées par moins de 10 % des séquences, il reste 2 886 colonnes en ADN-A et une seule séquence réellement lacunaire sur 974 (ADN-B : 3 292 colonnes, 27 séquences sur 601).

Deux alignements servent ensuite. Un jeu de 120 génomes par segment pour les détecteurs natifs, et un jeu réduit de 45 génomes pour la comparaison à sept méthodes, le coût des méthodes Python d'OpenRDP interdisant de les lancer sur 120. Dans les deux cas, les accessions de référence et des représentants de chaque groupe viral cité dans la table de vérité sont forcés dans l'alignement : sans cela l'écrémage par points les plus éloignés supprimait précisément les séquences sur lesquelles porte la notation.

## La table de vérité, et pourquoi elle n'est pas indépendante

25 événements ont été recensés : 16 sur l'ADN-A, 9 sur l'ADN-B. L'essentiel vient du tableau 1 de De Bruyn et al. 2012 (*BMC Evolutionary Biology* 12:228), complété par l'événement recombinant d'ACMBFV décrit par Tiendrébéogo et al. 2012 (*Virology Journal* 9:67). Parmi ces 25, 19 sont notables — 15 en ADN-A, 4 en ADN-B ; les autres désignent un recombinant trop large (« EACMV-like » parent de lui-même) ou un virus de tomate absent d'un jeu manioc.

Les coordonnées publiées sont exprimées dans le repère d'une accession, pas d'un alignement. Pour l'ADN-B, le texte de De Bruyn le dit explicitement pour l'événement B02 : région intergénique entre BV1 et BC1, positions 1124 à 1454 de JF909200. Ces bornes reproduisent exactement la ligne du tableau, ce qui fixe le repère pour tout le tableau ADN-B. Le même raisonnement vaut pour l'ADN-A avec JF909063, dont la fenêtre CP centrale (532–1014) encadre de peu les bornes de l'événement A01 (545–1008).

La projection vers les colonnes de nos alignements passe par la rotation appliquée à chaque génome, récupérée en recherchant la séquence orientée dans la concaténation de la séquence brute et de son complément inverse plutôt qu'en faisant confiance à un journal. Les trois accessions de référence donnent indépendamment le même décalage de 7 nucléotides, ce qui confirme que le recadrage est systématique. Contrôle biologique : l'ORF AV1 de JF909063 orienté occupe 342–1115 (257 acides aminés) et l'événement A01, projeté en 552–1015 dans ce repère, tombe entièrement à l'intérieur. C'est bien la région de la capside, celle que Zhou et al. avaient identifiée en 1997 sur des arguments d'identité de séquence.

Reste la limite de fond, qu'il faudra écrire noir sur blanc dans le manuscrit. Les 24 événements issus de De Bruyn ont été appelés avec RDP3, c'est-à-dire avec les sept méthodes que nous testons ici, et les bornes ont été « vérifiées et ajustées manuellement ». Mesurer 3SEQ ou BootScan contre ce tableau ne mesure pas leur exactitude : cela mesure leur capacité à reproduire un consensus auquel ils ont eux-mêmes contribué. Un échec est donc informatif ; une réussite ne l'est qu'à moitié.

## Métriques, et le témoin qu'il a fallu ajouter

Les définitions reprennent celles du banc de simulation : un événement est *signalé* si une méthode désigne un membre du groupe recombinant avec p < 0,05, et *localisé* si l'intervalle rapporté atteint un indice de Jaccard circulaire de 0,3 avec l'intervalle publié.

Appliquées telles quelles, ces définitions donnent un résultat qui aurait été trompeur. Toutes les méthodes signalent 100 % des événements, sur les deux segments, dans toutes les configurations — 3SEQ désigne les 120 génomes sur 120 comme recombinants en ADN-A. Et RDP localise 15 événements sur 15 : en émettant 26 219 intervalles sur 45 séquences. À ce volume, recouvrir une cible n'est plus une performance.

D'où un témoin de placement aléatoire, décidé après ce constat mais appliqué uniformément. Pour chaque événement et chaque méthode, on conserve le nombre d'appels portant sur le groupe recombinant et la longueur de chacun, et on retire au hasard leur position sur l'alignement circulaire, 200 fois. On en tire le nombre d'événements qu'une méthode de même volume localiserait par hasard, et une valeur-p empirique sur le meilleur recouvrement observé.

## État des détecteurs

![Temps de calcul]({{artifact:art_a309b759-1bc1-4960-8451-126ba0279054}})

Trois des cinq méthodes Python d'OpenRDP ne sont pas utilisables sur ce jeu.

MaxChi s'arrête sur `IndexError: index 2886 is out of bounds` à `openrdp/maxchi.py:264`. La fenêtre glissante s'étend d'un cran au-delà du dernier site sans contrôle de borne. La rustine fournie (`validation/openrdp_maxchi_bounds.patch`, contre le commit b42d9e5) lève ce point, et découvre immédiatement un second échec : `chi2_contingency` refuse une table dont un effectif attendu est nul.

Chimaera s'arrête sur `ZeroDivisionError` à `chimaera.py:238`, sur un `len(comp_seq)` nul. La ligne fautive mérite un second regard : le seuil s'y écrit `0.05/len(comp_seq) - win_size`, soit une probabilité moins une longueur de fenêtre. Le seuil de Bonferroni est mal formé indépendamment du plantage.

SiScan n'a pas terminé en 30 minutes sur 30 génomes, ni en 2 heures sur 45. Les deux lancements ont été interrompus.

RDP et BootScan aboutissent, en 537 et 1 995 secondes sur l'ADN-A réduit. Rapportés à 3SEQ sur le même alignement, cela fait 384 fois et 1 450 fois. La colonne `Pvalue` que RDP écrit atteint 13 704 et n'est donc pas une probabilité ; le code d'OpenRDP le reconnaît lui-même dans un commentaire `TODO` à `rdp.py:36`. La seule parallélisation prévue par OpenRDP passe par MPI, qui ne s'initialise pas dans l'environnement conteneurisé utilisé ici (découverte de topologie hwloc indisponible).

OpenRDP ne peut donc pas remplacer RDP4 en l'état. Deux méthodes sur sept, dont une qui ne rapporte pas de p-value exploitable, ce n'est pas un substitut.

## Résultats

![Localisation contre témoin]({{artifact:art_d01d1d75-75f5-47a5-a7ee-d5e9c961be54}})

Aucune méthode prise isolément ne localise les événements publiés mieux qu'un tirage au hasard de même volume. En ADN-A, 3SEQ en localise 3 contre 4,3 attendus par hasard, GENECONV 6 contre 10,0, BootScan 14 contre 15,0, RDP 15 contre 15,0. La valeur-p empirique médiane vaut 0,95 pour 3SEQ et 1,00 pour GENECONV, BootScan et RDP. Sur les 19 événements notables croisés avec les quatre méthodes, soit 76 combinaisons, quatre seulement passent sous 0,05 — toutes sur l'ADN-A.

Les parents sont un échec net. Sur l'ADN-A, cinq événements avaient leurs deux groupes parentaux présents dans l'alignement et une localisation par 3SEQ : A01 (EACMV-UG), A04, A06, A09 et A11. Dans aucun des cinq 3SEQ ne désigne un parent dans le groupe mineur et un dans le groupe majeur publiés. Zéro sur cinq.

![Carte de l'ADN-A]({{artifact:art_cd876023-78b6-4e7d-9b6c-cb788723ac0e}})

La carte explique pourquoi. Les points de rupture appelés couvrent la totalité du génome pour toutes les méthodes, avec des pics autour de la fin de l'AV1 et dans l'AC1, c'est-à-dire là où la littérature situe les points chauds — mais sans creux ailleurs. Il n'y a pas de signal à extraire d'un tel profil par un simple critère de recouvrement.

![Volume contre coïncidence]({{artifact:art_41690314-5a1a-4a30-b31c-bbd06c68ece6}})

C'est ici que la concordance paye. Exiger qu'au moins trois détecteurs s'accordent sur le même événement ramène 71 343 appels à 166 sur l'ADN-A et 39 312 à 94 sur l'ADN-B. La proportion d'appels qui coïncident avec un événement publié passe à 52,4 % et 46,8 %, contre 45,0 % et 9,0 % pour BootScan, 35,6 % et 6,0 % pour RDP. Soit 158 à 268 fois moins d'appels que les méthodes les plus prolixes, pour un taux de coïncidence supérieur sur les deux segments.

La concordance ne fait toujours pas mieux que le témoin aléatoire sur le critère de localisation : 9 événements sur 15 contre 9,6 attendus en ADN-A. Le seul endroit où l'observation dépasse le témoin est l'ADN-B à trois outils, 4 contre 2,8, ce qui ne veut pas dire grand-chose avec quatre événements.

## Ce que cela dit pour le manuscrit

Le résultat central n'est pas celui qui était attendu, et il vaut mieux qu'un résultat attendu. Sur un jeu réel de begomovirus du manioc, les sept méthodes classiques émettent des appels dont le recouvrement avec les événements publiés s'explique entièrement par leur volume. Ce n'est pas un argument contre les méthodes : c'est un argument contre leur lecture naïve, celle qui consiste à prendre une sortie de RDP et à en lire les intervalles comme des faits. La valeur de vipan se formule alors sans exagération : réduire un bruit de 26 000 à 44 000 appels à quelques centaines, puis soumettre ce qui reste à une validation phylogénétique.

Deux chiffres sont directement citables en introduction ou en discussion : 100 % des événements signalés par toutes les méthodes sur les deux segments, et zéro paire de parents correcte sur cinq pour 3SEQ. Le second est le plus parlant pour un virologue, parce que c'est le parent qui porte l'interprétation biologique.

Le banc de simulation reste indispensable au papier. Il est le seul endroit où la vérité est exacte, et c'est lui qui montre que ces méthodes fonctionnent quand elles sont nourries d'un recombinant à points de rupture connus. La section manioc vient après, et dit ce qui se passe quand on les lâche sur un jeu réel.

## Limites

RDP4 n'a pas pu être comparé : c'est une application Windows à interface graphique, hors de portée d'une exécution automatisée. La comparaison porte donc sur OpenRDP, dont la maturité est elle-même un résultat de l'étude.

La table de vérité n'est pas exhaustive. Elle recense ce qu'une étude de 2012 a choisi de rapporter, pas tous les événements présents dans un jeu qui compte aujourd'hui 974 ADN-A. Un appel sans correspondance dans la table n'est donc pas un faux positif, et aucune précision n'est calculée ici — seulement un taux de coïncidence, qui est autre chose.

L'ADN-B ne pèse que quatre événements notables, et leurs intervalles sont énormes : B01 couvre 43 % de l'alignement. N'importe quel appel tombant dans cette fenêtre compte comme une coïncidence. Les chiffres ADN-B sont à lire avec cette réserve.

Les 45 génomes du jeu réduit ne sont pas un échantillon neutre : ils contiennent par construction les recombinants et parents de référence, ce qui favorise plutôt la détection.

Dernier point, et c'est une contrainte de reproductibilité plus qu'une limite scientifique. MAFFT `--auto` ne rend pas le même alignement d'un passage à l'autre sur ce type de jeu : sur les 465 génomes sweepovirus du même projet, le même FASTA donnait 3 982 ou 3 969 colonnes selon le lancement. Toute la projection de la vérité-terrain repose sur des numéros de colonnes. Les quatre alignements utilisés ici sont donc déposés tels quels (`cassava_segA_alignment_final.fasta` et ses trois homologues, sommes MD5 dans le journal de session) et c'est sur eux qu'il faut rejouer les calculs, pas sur un alignement refait.

## Refaire les calculs

```bash
cd ~/vipan_transfer
./cassava/run_cassava.sh A          # et B
tools/venv/bin/python validation/make_keep_list.py --manifest work/cassava/raw/manifest.tsv \
  --truth validation/truth_cassava.tsv --oriented work/cassava/segA/oriented.fasta \
  --segment A --out-dir validation
tools/venv/bin/python validation/project_truth.py --truth validation/truth_cassava.tsv --segment A \
  --raw work/cassava/raw/raw.fasta --oriented work/cassava/segA/oriented.fasta \
  --aligned work/cassava/segA/aligned.fasta \
  --columns work/cassava/segA/alignment_final.fasta.columns.tsv \
  --out validation/truth_aligned_segA.tsv
tools/venv/bin/python validation/concordance.py --events <events.tsv ...> \
  --aln work/cassava/segA/alignment_reduced.fasta --out work/cassava/segA/concordance/events.tsv
tools/venv/bin/python validation/score_cassava.py --truth validation/truth_cassava.tsv \
  --truth-aligned validation/truth_aligned_segA.tsv --groups validation/groups_segA.tsv \
  --aln work/cassava/segA/alignment_reduced.fasta --segment A \
  --events <events.tsv ...> --out-dir validation/scores_reduced
```

Les scripts ajoutés sont dans `validation/` : `truth_cassava.tsv` (table de vérité, une ligne par événement avec sa source), `make_keep_list.py`, `project_truth.py`, `openrdp_to_events.py`, `concordance.py`, `score_cassava.py`, `timings.tsv` et `openrdp_maxchi_bounds.patch`. Le filtre corrigé est `cassava/postalign.py`.
