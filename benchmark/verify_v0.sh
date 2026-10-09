#!/usr/bin/env bash
# Compares the freshly generated truth tables and scores with the ones from the first run (expected_v0/).
set -u
cd "$(dirname "$0")"
ok=1
for pair in "real:work/data" "synthetic:work/data_syn"; do
  e=${pair%%:*}; w=${pair##*:}
  cmp -s expected_v0/$e/truth.tsv $w/truth.tsv && echo "truth $e : identique" || { echo "truth $e : DIFFERENT"; ok=0; }
done
for pair in "real:work/out" "synthetic:work/out_syn"; do
  e=${pair%%:*}; w=${pair##*:}
  cmp -s expected_v0/$e/summary.txt $w/summary.txt && echo "summary $e : identique" || { echo "summary $e : DIFFERENT (voir diff)"; ok=0; }
done
exit $((1-ok))
