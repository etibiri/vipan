#!/usr/bin/env bash
# Prepare tools: OpenRDP at a pinned commit (it ships the 3SEQ, GENECONV binaries and the 3SEQ p-value table).
# Needs: git, python3 (3.9 to 3.12 recommended), internet access to github.com. x86-64 Linux only.
set -euo pipefail
cd "$(dirname "$0")"
[ "$(uname -m)" = "x86_64" ] || { echo "Les binaires 3SEQ/GENECONV fournis par OpenRDP sont x86-64 uniquement."; exit 1; }
mkdir -p tools work
if [ ! -d tools/OpenRDP ]; then
  git clone https://github.com/PoonLab/OpenRDP.git tools/OpenRDP
fi
git -C tools/OpenRDP checkout b42d9e5fb84c14e817e08aab2668aeb2c6a94f83
chmod +x tools/OpenRDP/openrdp/bin/3Seq/3seq.Unix tools/OpenRDP/openrdp/bin/GENECONV/geneconv.Unix
python3 -m venv tools/venv
tools/venv/bin/pip install --quiet -r requirements.txt
# sanitised copy of the alignment (OpenRDP and the native tools reject IUPAC codes other than N)
mkdir -p work/data
tools/venv/bin/python scripts/sanitize.py data/final_alignment_316.fasta work/data/full316_sanitised.fa
echo "setup ok"
echo "Optionnel (methodes Python d'OpenRDP : RDP, BootScan, SiScan, Chimaera, MaxChi) : python 3.9-3.12 requis, puis"
echo "  tools/venv/bin/pip install -r requirements-openrdp-python.txt && tools/venv/bin/pip install --no-deps tools/OpenRDP"
