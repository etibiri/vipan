#!/usr/bin/env python3
"""Convert the CSV written by OpenRDP's Python methods into the common events.tsv format.

OpenRDP writes Method,Start,End,Recombinant,Parent1,Parent2,Pvalue and applies its own per-method
significance filter before writing, so no threshold is re-applied here.

One quirk has to be handled rather than silently propagated: the Rdp method reports values far above 1 in
the Pvalue column (13704 on our data), which cannot be a probability. OpenRDP's own source carries a TODO
on exactly this point ("no innate max_pvalue, output has given values greater than 31 before",
openrdp/rdp.py line 36). For Rdp rows the reported number is therefore kept in `note` and `p_corrected` is
set to 0, meaning "call retained, significance as judged by OpenRDP itself". The other four methods emit
probabilities below 1 and are passed through unchanged.
"""
import argparse
import csv
import os

TOOL = {"rdp": "RDP", "bootscan": "BootScan", "maxchi": "MaxChi", "siscan": "SiScan", "chimaera": "Chimaera"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows, skipped = [], {}
    for path in a.csv:
        method = os.path.splitext(os.path.basename(path))[0]
        tool = TOOL.get(method, method)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            skipped[tool] = "no output"
            continue
        n = 0
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if not r.get("Recombinant"):
                    continue
                val = float(r["Pvalue"])
                p = 0.0 if val > 1 else val
                note = f"reported={r['Pvalue']}" if val > 1 else ""
                rows.append((tool, r["Recombinant"], r["Parent1"], r["Parent2"],
                             int(float(r["Start"])), int(float(r["End"])), p, note))
                n += 1
        skipped[tool] = f"{n} calls"
    with open(a.out, "w") as f:
        f.write("tool\trecombinant\tparent_a\tparent_b\tstart\tend\tp_corrected\tnote\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"{a.out}: " + ", ".join(f"{k} {v}" for k, v in sorted(skipped.items())))


if __name__ == "__main__":
    main()
