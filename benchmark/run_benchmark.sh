#!/usr/bin/env bash
# Reproduces benchmark v0: simulated recombinants (real-alignment and synthetic backgrounds),
# 3SEQ and GENECONV, scoring. WORKERS = parallel jobs (default: all cores).
set -euo pipefail
cd "$(dirname "$0")"
export VIPAN_ROOT="$PWD"
PY=tools/venv/bin/python
mkdir -p work/data work/data_syn work/out work/out_syn
$PY scripts/sim_make.py
$PY scripts/sim_synth.py
$PY scripts/run_detect.py                      > work/out/run.log
SIM_SET=_syn $PY scripts/run_detect.py         > work/out_syn/run.log
$PY scripts/analyse.py                         > work/out/summary_stdout.txt
SIM_SET=_syn $PY scripts/analyse.py            > work/out_syn/summary_stdout.txt
echo "Resultats : work/out/summary.txt (fond reel) et work/out_syn/summary.txt (fond synthetique)"
