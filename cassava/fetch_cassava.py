#!/usr/bin/env python3
"""Download cassava mosaic begomovirus genomes from NCBI nuccore (E-utilities, standard library only).

Writes, in --out:
  raw.fasta       every record kept after the length and "complete" filters
  manifest.tsv    accession, organism, isolate, country, length, segment call (A, B or unassigned), title
Set NCBI_EMAIL (and NCBI_API_KEY if you have one) in the environment: NCBI asks for it.
"""
import argparse, json, os, re, sys, time
import urllib.error, urllib.parse, urllib.request

EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SPECIES = [
    "African cassava mosaic virus",
    "East African cassava mosaic virus",
    "East African cassava mosaic Cameroon virus",
    "East African cassava mosaic Kenya virus",
    "East African cassava mosaic Malawi virus",
    "East African cassava mosaic Zanzibar virus",
    "African cassava mosaic Burkina Faso virus",
    "South African cassava mosaic virus",
    "Cassava mosaic Madagascar virus",
    "Indian cassava mosaic virus",
    "Sri Lankan cassava mosaic virus",
]
RE_A = re.compile(r"\bDNA[- ]?A\b|\bsegment A\b|\bcomponent A\b", re.I)
RE_B = re.compile(r"\bDNA[- ]?B\b|\bsegment B\b|\bcomponent B\b", re.I)
RE_COMPLETE = re.compile(r"complete (genome|sequence|nucleotide)", re.I)


def call(endpoint, params, retries=6):
    params = dict(params)
    if os.environ.get("NCBI_EMAIL"):
        params["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    params["tool"] = "vipan-cassava-fetch"
    data = urllib.parse.urlencode(params).encode()
    for k in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(f"{EUT}/{endpoint}", data=data), timeout=60) as r:
                return r.read().decode()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            wait = 2 ** k
            print(f"  {endpoint}: {e}; retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise SystemExit(f"{endpoint} failed after {retries} attempts")


def esearch(term, retmax=5000):
    out = json.loads(call("esearch.fcgi", dict(db="nuccore", term=term, retmax=retmax, retmode="json")))
    return out["esearchresult"]["idlist"]


def parse_summary(js):
    """esummary JSON -> list of dicts."""
    res = js["result"]
    rows = []
    for uid in res.get("uids", []):
        d = res[uid]
        subtype = (d.get("subtype") or "").split("|")
        subname = (d.get("subname") or "").split("|")
        sub = dict(zip(subtype, subname))
        rows.append(dict(uid=uid, accession=d.get("accessionversion") or d.get("caption", ""),
                         organism=d.get("organism", ""), title=d.get("title", ""),
                         length=int(d.get("slen") or 0), isolate=sub.get("isolate", ""),
                         country=sub.get("country", ""), collection_date=sub.get("collection_date", "")))
    return rows


def segment_call(title):
    a, b = bool(RE_A.search(title)), bool(RE_B.search(title))
    if a and not b:
        return "A"
    if b and not a:
        return "B"
    return "unassigned"


def esummary(ids):
    rows = []
    for i in range(0, len(ids), 100):
        rows += parse_summary(json.loads(call("esummary.fcgi", dict(db="nuccore", id=",".join(ids[i:i + 100]), retmode="json"))))
        time.sleep(0.4)
    return rows


def efetch_fasta(accessions):
    chunks = []
    for i in range(0, len(accessions), 100):
        chunks.append(call("efetch.fcgi", dict(db="nuccore", id=",".join(accessions[i:i + 100]), rettype="fasta", retmode="text")))
        time.sleep(0.4)
    return "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="cassava_data")
    ap.add_argument("--min-len", type=int, default=2400)
    ap.add_argument("--max-len", type=int, default=3000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    seen, rows = set(), []
    for sp in SPECIES:
        ids = esearch(f'"{sp}"[Organism] AND {a.min_len}:{a.max_len}[SLEN]')
        time.sleep(0.4)
        new = [i for i in ids if i not in seen]
        seen.update(new)
        print(f"{sp}: {len(ids)} records ({len(new)} new)")
        rows += esummary(new) if new else []
    kept = [r for r in rows if RE_COMPLETE.search(r["title"])]
    print(f"{len(rows)} records, {len(kept)} with 'complete genome/sequence' in the title")
    for r in kept:
        r["segment"] = segment_call(r["title"])
    fasta = efetch_fasta([r["accession"] for r in kept])
    open(f"{a.out}/raw.fasta", "w").write(fasta)
    cols = ["accession", "organism", "isolate", "country", "collection_date", "length", "segment", "title"]
    with open(f"{a.out}/manifest.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in kept:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    n = {s: sum(r["segment"] == s for r in kept) for s in ("A", "B", "unassigned")}
    print("segments:", n, "->", a.out)


if __name__ == "__main__":
    main()
