#!/usr/bin/env python3
"""Score cassava recombination calls against the published truth table.

Definitions, fixed before looking at any result and deliberately parallel to those of the simulation bench
(scripts/analyse.py), with one change: the truth is a virus group, not a single simulated sequence, because
the published events are attributed to species or strain groups rather than to individual accessions.

  flagged     : at least one significant call (p < 0.05) names a member of the event's recombinant group as
                recombinant. 3SEQ names a child; GENECONV only names a pair, so either member counts.
  localised   : flagged, and the reported interval reaches a circular Jaccard index of at least 0.3 with the
                published interval.
  parents ok  : among the localised 3SEQ calls, at least one reports two parents that fall in the published
                minor and major parent groups, in either order. Scored only when both parent groups are
                named and present in the alignment.
  bp error    : for the best-Jaccard localised call, the distance in alignment columns between reported and
                published breakpoints.

What is NOT measured here: precision. The truth table lists the events a 2012 study chose to report, not
every event present in this dataset, so a call that matches nothing in the table is not a false positive. It
is counted and reported as an unmatched call, which is a different quantity.
"""
import argparse
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_keep_list import EACMV_LIKE, SPECIFIC_ACCESSIONS  # noqa: E402

ALPHA, JAC = 0.05, 0.3


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def aln_length(path):
    n = 0
    for line in open(path):
        if line.startswith(">"):
            if n:
                return n
        else:
            n += len(line.strip())
    return n


def arc(interval, L):
    """Interval (begin, end) in 1-based columns -> (start, length) on the circle, wrapping if begin > end."""
    b, e = interval
    return b, (e - b) % L + 1


def jaccard(i1, i2, L):
    """Jaccard index of two arcs on a circular alignment, computed arithmetically.

    The set-of-columns version is equivalent but far too slow for the permutation control, which evaluates
    this millions of times.
    """
    a0, la = arc(i1, L)
    b0, lb = arc(i2, L)
    inter = 0
    for k in (-L, 0, L):
        lo, hi = max(a0, b0 + k), min(a0 + la, b0 + k + lb)
        if hi > lo:
            inter += hi - lo
    union = la + lb - inter
    return inter / union if union > 0 else 0.0


def circ_dist(x, y, L):
    d = abs(x - y)
    return min(d, L - d)


def members(group, acc2grp):
    """Accessions of a truth-table group present in the alignment."""
    if group in SPECIFIC_ACCESSIONS:
        return {a for a in SPECIFIC_ACCESSIONS[group] if a in acc2grp}
    if group.startswith("EACMV-[") or group.startswith("ACMBFV-["):
        return {group} & set(acc2grp)
    out = set()
    for g in group.split("|"):
        g = g.strip()
        if g == "EACMV-like":
            out |= {a for a, gg in acc2grp.items() if gg in EACMV_LIKE}
        elif g and g != "Unknown":
            out |= {a for a, gg in acc2grp.items() if gg == g}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--truth-aligned", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--segment", required=True)
    ap.add_argument("--events", required=True, nargs="+",
                    help="events.tsv files; all tools found inside are scored separately")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--permutations", type=int, default=200,
                    help="random-placement control replicates per event; 0 disables")
    ap.add_argument("--seed", type=int, default=20261009)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    rng_np = np.random.default_rng(a.seed)
    L = aln_length(a.aln)
    in_aln = {l[1:].split()[0] for l in open(a.aln) if l.startswith(">")}
    acc2grp = {r["accession"]: r["group"] for r in read_tsv(a.groups) if r["accession"] in in_aln}

    truth = {r["event_id"]: r for r in read_tsv(a.truth) if r["segment"] == a.segment}
    for r in read_tsv(a.truth_aligned):
        if r["event_id"] in truth:
            truth[r["event_id"]].update({"begin_col": r["begin_col"], "end_col": r["end_col"]})

    calls = []
    for path in a.events:
        for r in read_tsv(path):
            if float(r["p_corrected"]) >= ALPHA:
                continue
            calls.append({"tool": r["tool"], "recombinant": r["recombinant"],
                          "parents": (r["parent_a"], r["parent_b"]),
                          "interval": (int(r["start"]), int(r["end"])),
                          "p": float(r["p_corrected"])})
    tools = sorted({c["tool"] for c in calls})

    rows, summary, matched_calls = [], [], set()
    for tool in tools:
        tcalls = [c for c in calls if c["tool"] == tool]
        names_recomb = tool == "3SEQ"
        per_event = []
        for eid in sorted(truth):
            t = truth[eid]
            if t["scorable"] != "yes" or t.get("begin_col", "NA") == "NA":
                continue
            grp = members(t["recombinant_group"], acc2grp)
            ti = (int(t["begin_col"]), int(t["end_col"]))
            if names_recomb:
                hits = [c for c in tcalls if c["recombinant"] in grp]
            else:
                hits = [c for c in tcalls if set(c["parents"]) & grp]
            best, bestj = None, 0.0
            for c in hits:
                j = jaccard(ti, c["interval"], L)
                if j > bestj:
                    best, bestj = c, j
            minor, major = members(t["minor_parent_group"], acc2grp), members(t["major_parent_group"], acc2grp)
            parents_ok = ""
            if names_recomb and minor and major:
                loc = [c for c in hits if jaccard(ti, c["interval"], L) >= JAC]
                parents_ok = str(any({c["parents"][0], c["parents"][1]} & minor and
                                     {c["parents"][0], c["parents"][1]} & major for c in loc))
            localised = bestj >= JAC
            # Random-placement control: keep this tool's calls on this group, keep their lengths, draw new
            # start positions uniformly on the circular alignment. A method that emits tens of thousands of
            # intervals localises any target by volume alone; this says how much of the observed
            # localisation that accounts for.
            p_chance, p_emp = "", ""
            if hits and a.permutations:
                lens = np.array([arc(c["interval"], L)[1] for c in hits], dtype=np.int64)
                t0, lt = arc(ti, L)
                S = rng_np.integers(1, L + 1, size=(a.permutations, lens.size))
                inter = np.zeros_like(S, dtype=np.int64)
                for k in (-L, 0, L):
                    lo = np.maximum(t0, S + k)
                    hi = np.minimum(t0 + lt, S + k + lens)
                    inter += np.maximum(0, hi - lo)
                jacs = inter / (lt + lens - inter)
                best_perm = jacs.max(axis=1)
                p_chance = round(float((best_perm >= JAC).mean()), 3)
                p_emp = round(float((best_perm >= bestj).mean()), 3)
            if localised:
                for c in hits:
                    if jaccard(ti, c["interval"], L) >= JAC:
                        matched_calls.add((c["tool"], c["recombinant"], c["parents"], c["interval"]))
            rows.append({
                "segment": a.segment, "tool": tool, "event_id": eid,
                "recombinant_group": t["recombinant_group"], "group_size": len(grp),
                "truth_begin": ti[0], "truth_end": ti[1],
                "n_calls_on_group": len(hits), "flagged": bool(hits), "localised": localised,
                "best_jaccard": round(bestj, 3),
                "call_begin": best["interval"][0] if best else "", "call_end": best["interval"][1] if best else "",
                "err_begin": circ_dist(best["interval"][0], ti[0], L) if (best and localised) else "",
                "err_end": circ_dist(best["interval"][1], ti[1], L) if (best and localised) else "",
                "parents_ok": parents_ok, "best_p": best["p"] if best else "", "localised_by_chance": p_chance, "p_empirical": p_emp,
            })
            per_event.append(rows[-1])
        n = len(per_event)
        fl = sum(r["flagged"] for r in per_event)
        lo = sum(r["localised"] for r in per_event)
        errs_b = [r["err_begin"] for r in per_event if r["err_begin"] != ""]
        errs_e = [r["err_end"] for r in per_event if r["err_end"] != ""]
        pa = [r["parents_ok"] for r in per_event if r["parents_ok"] != ""]
        chance = [r["localised_by_chance"] for r in per_event if r["localised_by_chance"] != ""]
        pemp = [r["p_empirical"] for r in per_event if r["p_empirical"] != ""]
        unmatched = [c for c in tcalls
                     if (c["tool"], c["recombinant"], c["parents"], c["interval"]) not in matched_calls]
        summary.append({
            "segment": a.segment, "tool": tool, "n_events": n,
            "flagged": fl, "flagged_pct": round(100 * fl / n, 1) if n else "",
            "localised": lo, "localised_pct": round(100 * lo / n, 1) if n else "",
            "localised_expected_by_chance": round(sum(chance), 2) if chance else "",
            "events_better_than_chance": sum(x <= 0.05 for x in pemp),
            "median_p_empirical": round(sorted(pemp)[len(pemp) // 2], 3) if pemp else "",
            "parents_ok": sum(x == "True" for x in pa), "parents_scored": len(pa),
            "median_err_begin": sorted(errs_b)[len(errs_b) // 2] if errs_b else "",
            "median_err_end": sorted(errs_e)[len(errs_e) // 2] if errs_e else "",
            "n_calls": len(tcalls), "n_calls_unmatched": len(unmatched),
        })

    os.makedirs(a.out_dir, exist_ok=True)
    def write(path, data, cols):
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in data:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")

    write(os.path.join(a.out_dir, f"per_event_seg{a.segment}.tsv"), rows, list(rows[0]) if rows else ["segment"])
    write(os.path.join(a.out_dir, f"summary_seg{a.segment}.tsv"), summary, list(summary[0]))
    for s in summary:
        print(f"seg{a.segment} {s['tool']:<10} {s['flagged']}/{s['n_events']} flagged, "
              f"{s['localised']}/{s['n_events']} localised "
              f"(chance {s['localised_expected_by_chance']}), "
              f"{s['events_better_than_chance']}/{s['n_events']} beat the control, {s['n_calls']} calls "
              f"({s['n_calls_unmatched']} unmatched)")


if __name__ == "__main__":
    main()
