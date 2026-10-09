#!/usr/bin/env bash
# Cassava mosaic begomoviruses: download -> orient -> MAFFT -> clean/thin -> native 3SEQ + GENECONV -> events.tsv
# Usage: cassava/run_cassava.sh [A|B]   (default A: DNA-A; the two segments must be analysed separately)
# Environment: NCBI_EMAIL (recommended), THREADS (default: all cores), MAX_SEQS (default 120, 0 = no cap),
#              MIN_LEN / MAX_LEN (default 2500 / 2900, accepted genome length)
set -euo pipefail
cd "$(dirname "$0")/.."
export VIPAN_ROOT="$PWD"
SEG="${1:-A}"
PY=tools/venv/bin/python
D="work/cassava/seg$SEG"
mkdir -p "$D"
command -v mafft >/dev/null || { echo "mafft manquant : sudo apt install mafft"; exit 1; }
if [ ! -s work/cassava/raw/raw.fasta ]; then
  $PY cassava/fetch_cassava.py --out work/cassava/raw
fi
$PY cassava/prepare_cassava.py --raw work/cassava/raw/raw.fasta --manifest work/cassava/raw/manifest.tsv --segment "$SEG" --min-len "${MIN_LEN:-2500}" --max-len "${MAX_LEN:-2900}" --out "$D"
[ -s "$D/oriented.fasta" ] || { echo "aucune sequence retenue, voir $D/prepare_log.tsv"; exit 1; }
mafft --auto --thread "${THREADS:-$(nproc)}" --inputorder "$D/oriented.fasta" > "$D/aligned.fasta" 2> "$D/mafft.log"
$PY cassava/postalign.py --aln "$D/aligned.fasta" --out "$D/alignment_final.fasta" --log "$D/thinning_log.tsv" --max-seqs "${MAX_SEQS:-120}"
$PY cassava/detect_native.py --aln "$D/alignment_final.fasta" --out "$D/detect"
echo "Resultats : $D/detect/events.tsv ; manifeste : work/cassava/raw/manifest.tsv"
