"""Replace IUPAC ambiguity codes other than N with N, as an adapter would before OpenRDP."""
import re
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
changed = 0
out = []
for line in src.read_text().splitlines():
    if line.startswith(">"):
        out.append(line)
        continue
    cleaned, n = re.subn(r"[^ACGTNacgtn-]", "N", line)
    changed += n
    out.append(cleaned)
dst.write_text("\n".join(out) + "\n")
print(src.name, "->", dst.name, "| ambiguity characters replaced:", changed)
