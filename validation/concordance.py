#!/usr/bin/env python3
"""Build event-level concordance call sets from several single-tool events.tsv files.

This is the step vipan exists for. Each detector on its own emits calls whose overlap with the published
events is no better than randomly placed intervals of the same sizes; the question is whether requiring
several detectors to agree on the same event concentrates the calls onto the real ones.

Clustering rule, fixed in advance:
  two calls belong to the same event if they implicate the same sequence and their intervals reach a
  circular Jaccard index of at least --min-jaccard (0.5 by default). "Implicate the same sequence" means
  the same named recombinant for tools that name one, and shared membership of the reported pair for
  GENECONV, which does not name a recipient.
  A cluster's support is its number of distinct tools. Its representative interval is the medoid call, the
  one with the highest summed Jaccard to the rest of the cluster.

Output is an events.tsv with tool = CONCORDANCE<k>, one row per cluster supported by at least k tools, so
it can be scored by exactly the same code and the same random-placement control as a single detector.
"""
import argparse
import os


def arc(interval, L):
    b, e = interval
    return b, (e - b) % L + 1


def jaccard(i1, i2, L):
    a0, la = arc(i1, L)
    b0, lb = arc(i2, L)
    inter = 0
    for k in (-L, 0, L):
        lo, hi = max(a0, b0 + k), min(a0 + la, b0 + k + lb)
        if hi > lo:
            inter += hi - lo
    union = la + lb - inter
    return inter / union if union > 0 else 0.0


def aln_length(path):
    n = 0
    for line in open(path):
        if line.startswith(">"):
            if n:
                return n
        else:
            n += len(line.strip())
    return n


def read_events(paths):
    out = []
    for p in paths:
        with open(p) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                out.append({"tool": r["tool"], "recombinant": r["recombinant"],
                            "parents": (r["parent_a"], r["parent_b"]),
                            "interval": (int(r["start"]), int(r["end"])),
                            "p": float(r["p_corrected"])})
    return out


def subject(call):
    """The sequences a call implicates as recombinant."""
    if call["recombinant"] not in ("NA", ""):
        return (call["recombinant"],)
    return call["parents"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True, nargs="+")
    ap.add_argument("--aln", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-jaccard", type=float, default=0.5)
    ap.add_argument("--support", type=int, nargs="+", default=[2, 3])
    a = ap.parse_args()

    L = aln_length(a.aln)
    calls = read_events(a.events)
    by_seq = {}
    for c in calls:
        for s in subject(c):
            by_seq.setdefault(s, []).append(c)

    clusters = []
    for seq, group in sorted(by_seq.items()):
        group = sorted(group, key=lambda c: (c["interval"], c["tool"]))
        assigned = [None] * len(group)
        for i, c in enumerate(group):
            if assigned[i] is not None:
                continue
            members = [i]
            assigned[i] = len(clusters)
            for j in range(i + 1, len(group)):
                if assigned[j] is None and jaccard(c["interval"], group[j]["interval"], L) >= a.min_jaccard:
                    assigned[j] = len(clusters)
                    members.append(j)
            clusters.append({"seq": seq, "calls": [group[k] for k in members]})

    rows = []
    for cl in clusters:
        tools = sorted({c["tool"] for c in cl["calls"]})
        if len(tools) < min(a.support):
            continue
        best, bestsum = None, -1.0
        for c in cl["calls"]:
            s = sum(jaccard(c["interval"], o["interval"], L) for o in cl["calls"])
            if s > bestsum:
                best, bestsum = c, s
        rows.append((cl["seq"], tools, best, min(c["p"] for c in cl["calls"])))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        f.write("tool\trecombinant\tparent_a\tparent_b\tstart\tend\tp_corrected\tnote\n")
        counts = {}
        for k in sorted(a.support):
            for seq, tools, best, p in rows:
                if len(tools) < k:
                    continue
                f.write("\t".join([f"CONCORDANCE{k}", seq, best["parents"][0], best["parents"][1],
                                   str(best["interval"][0]), str(best["interval"][1]), str(p),
                                   "+".join(tools)]) + "\n")
                counts[k] = counts.get(k, 0) + 1
    print(f"{a.out}: {len(calls)} calls -> {len(clusters)} clusters; " +
          ", ".join(f"support>={k}: {counts.get(k, 0)}" for k in sorted(a.support)))


if __name__ == "__main__":
    main()
