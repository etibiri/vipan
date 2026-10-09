"""Build seeded alignment subsets that always contain one positive-control event."""
import random
import sys
from pathlib import Path

src, outdir = Path(sys.argv[1]), Path(sys.argv[2])
must_keep = ["AJ132548.1", "DI183351.1", "MW574024.1"]

records, name, buf = {}, None, []
for line in src.read_text().splitlines():
    if line.startswith(">"):
        if name:
            records[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
records[name] = "".join(buf)

missing = [a for a in must_keep if a not in records]
print("records:", len(records), "| alignment length:", len(next(iter(records.values()))), "| missing controls:", missing)
if missing:
    sys.exit("control accessions not found; headers look like: " + ", ".join(list(records)[:3]))

others = sorted(set(records) - set(must_keep))
rng = random.Random(7)
rng.shuffle(others)
for n in (10, 20, 40, 80):
    chosen = must_keep + others[: n - len(must_keep)]
    with open(outdir / f"subset_{n}.fasta", "w") as handle:
        for acc in chosen:
            handle.write(f">{acc}\n{records[acc]}\n")
    print("wrote subset", n)
