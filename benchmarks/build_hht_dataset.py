#!/usr/bin/env python3
"""
HHT variant benchmark: text2gene2 vs LitVar2.

HHT genes: ENG (HHT1), ACVRL1 (HHT2), SMAD4 (JP-HHT), GDF2/BMP9 (HHT5)

Data sources (all bulk/local — no per-variant API calls for reference data):
  - ClinVar variant_summary.txt.gz + var_citations.txt  → variant list + ClinVar PMID counts
  - LitVar2 litvar2_variants.json.gz (streamed, HHT genes only) → LitVar2 PMID sets

text2gene2 is queried concurrently via its API (async, rate-limited).

Usage:
  python benchmarks/build_hht_dataset.py [--api http://loki.local] [--concurrency 8]
"""

import argparse
import asyncio
import csv
import gzip
import io
import json
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import httpx

HHT_GENES = {"ENG", "ACVRL1", "SMAD4", "GDF2"}

VARIANT_SUMMARY_URL  = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
VAR_CITATIONS_URL    = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"
LITVAR2_BULK_URL     = "https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz"

_HGVS_RE = re.compile(r"(N[MR]_\d+\.\d+)(?:\([^)]+\))?:(c\.[^\s(,]+)")


def parse_hgvs(name: str) -> str | None:
    m = _HGVS_RE.search(name)
    return f"{m.group(1)}:{m.group(2)}" if m else None


# ── Step 1: ClinVar citation counts ──────────────────────────────────────────

def fetch_citation_counts() -> dict[int, set[int]]:
    """AlleleID → set of PMIDs from ClinVar's explicit citations."""
    print("Step 1: ClinVar citation counts (var_citations.txt)...", flush=True)
    cites: dict[int, set[int]] = defaultdict(set)
    with urllib.request.urlopen(VAR_CITATIONS_URL) as resp:
        for line in io.TextIOWrapper(resp, encoding="utf-8", errors="replace"):
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6 and parts[4] == "PubMed":
                try:
                    cites[int(parts[0])].add(int(parts[5]))
                except ValueError:
                    pass
    print(f"  {len(cites):,} alleles have ≥1 ClinVar citation", flush=True)
    return cites


# ── Step 2: HHT variants from ClinVar ────────────────────────────────────────

def fetch_hht_variants(citation_map: dict[int, set[int]]) -> list[dict]:
    print("Step 2: HHT variants from variant_summary...", flush=True)
    variants, seen_hgvs = [], set()

    with urllib.request.urlopen(VARIANT_SUMMARY_URL) as resp:
        f = gzip.open(resp, "rt", encoding="utf-8", errors="replace")
        header = next(f).lstrip("#").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}

        def get(parts, col):
            i = idx.get(col)
            return parts[i].strip() if i is not None and i < len(parts) else ""

        for line in f:
            parts = line.rstrip("\n").split("\t")
            gene = get(parts, "GeneSymbol").split(";")[0].split(",")[0].strip()
            if gene not in HHT_GENES:
                continue
            if get(parts, "Assembly") != "GRCh38":
                continue
            if get(parts, "ClinSigSimple") != "1":   # P/LP only
                continue

            hgvs = parse_hgvs(get(parts, "Name"))
            if not hgvs or hgvs in seen_hgvs:
                continue
            seen_hgvs.add(hgvs)

            try:
                allele_id = int(get(parts, "AlleleID"))
            except ValueError:
                continue

            rsid = get(parts, "RS# (dbSNP)")
            variants.append({
                "hgvs":            hgvs,
                "gene":            gene,
                "disease":         get(parts, "PhenotypeList").split("|")[0],
                "clinsig":         get(parts, "ClinicalSignificance"),
                "allele_id":       allele_id,
                "variation_id":    get(parts, "VariationID"),
                "rsid":            rsid if rsid != "-1" else "",
                "review_status":   get(parts, "ReviewStatus"),
                "clinvar_pmids":   citation_map.get(allele_id, set()),
            })

    variants.sort(key=lambda v: (v["gene"], v["hgvs"]))
    print(f"  {len(variants)} P/LP HHT variants found across {HHT_GENES}", flush=True)
    return variants


# ── Step 3: LitVar2 PMID sets from bulk file ─────────────────────────────────

def build_litvar2_lookup(variants: list[dict]) -> dict[str, set[int]]:
    """
    Stream litvar2_variants.json.gz and build rsID → PMIDs + short-hgvs → PMIDs.

    LitVar2 all_hgvs format: "count|short_hgvs" e.g. "5|c.652C>T" or "6|p.R218W"
    — no transcript accession, just the change portion.
    We match our full HGVS (NM_xxx:c.yyy) by extracting the c./p. part.
    """
    print("Step 3: LitVar2 bulk data (streaming for HHT genes)...", flush=True)

    # Build short-HGVS index from our variants: "c.652C>T" → full hgvs
    short_to_full: dict[str, str] = {}
    for v in variants:
        if ":" in v["hgvs"]:
            short = v["hgvs"].split(":", 1)[1]   # "c.652C>T"
            short_to_full[short] = v["hgvs"]

    rsid_index: set[str] = {v["rsid"] for v in variants if v["rsid"]}

    lv2_by_full_hgvs: dict[str, set[int]] = defaultdict(set)
    lv2_by_rsid:      dict[str, set[int]] = defaultdict(set)

    n_records = 0
    with urllib.request.urlopen(LITVAR2_BULK_URL) as resp:
        f = gzip.open(resp, "rt", encoding="utf-8", errors="replace")
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue

            n_records += 1
            if n_records % 500_000 == 0:
                print(f"  ...{n_records/1_000_000:.1f}M LitVar2 records scanned", flush=True)

            gene_list = rec.get("gene") or []
            if isinstance(gene_list, str):
                gene_list = [gene_list]
            if not (set(gene_list) & HHT_GENES):
                continue

            pmids = set(int(p) for p in (rec.get("pmids") or []))
            if not pmids:
                continue

            # rsID match
            rsid = rec.get("rsid")
            if rsid and rsid in rsid_index:
                lv2_by_rsid[rsid].update(pmids)

            # all_hgvs match — format is "count|short_hgvs"
            for entry in (rec.get("all_hgvs") or []):
                short = entry.split("|", 1)[-1] if "|" in entry else entry
                if short in short_to_full:
                    lv2_by_full_hgvs[short_to_full[short]].update(pmids)

    print(f"  Scanned {n_records:,} records. "
          f"Matched {len(lv2_by_full_hgvs)} by HGVS, {len(lv2_by_rsid)} by rsID.", flush=True)

    # Merge per variant
    result: dict[str, set[int]] = {}
    for v in variants:
        pmids: set[int] = set()
        pmids.update(lv2_by_full_hgvs.get(v["hgvs"], set()))
        if v["rsid"]:
            pmids.update(lv2_by_rsid.get(v["rsid"], set()))
        result[v["hgvs"]] = pmids
    return result


# ── Step 4: text2gene2 concurrent queries ────────────────────────────────────

async def query_t2g2(hgvs: str, api_base: str,
                     client: httpx.AsyncClient,
                     sem: asyncio.Semaphore) -> set[int]:
    async with sem:
        try:
            r = await client.get(f"{api_base}/api/v2/hgvs2pmid/{hgvs}", timeout=90)
            return set(c["pmid"] for c in r.json().get("citations", []))
        except Exception as e:
            print(f"  t2g2 error {hgvs}: {e}", file=sys.stderr, flush=True)
            return set()


async def query_all_t2g2(variants: list[dict], api_base: str,
                          concurrency: int) -> dict[str, set[int]]:
    print(f"Step 4: querying text2gene2 ({concurrency} concurrent, {len(variants)} variants)...",
          flush=True)
    sem = asyncio.Semaphore(concurrency)
    limits = httpx.Limits(max_connections=concurrency, max_keepalive_connections=concurrency)
    async with httpx.AsyncClient(limits=limits) as client:
        tasks = {v["hgvs"]: query_t2g2(v["hgvs"], api_base, client, sem)
                 for v in variants}
        results = await asyncio.gather(*tasks.values())
    print(f"  Done.", flush=True)
    return dict(zip(tasks.keys(), results))


# ── Output ────────────────────────────────────────────────────────────────────

def write_tsv(rows: list[dict], path: Path) -> None:
    fields = ["hgvs", "gene", "disease", "clinsig",
              "n_clinvar", "n_t2g2", "n_lv2",
              "overlap", "only_t2g2", "only_lv2",
              "rsid", "review_status", "allele_id", "variation_id"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def print_summary(rows: list[dict]) -> None:
    print("\n" + "="*80)
    print(f"{'GENE':<10} {'HGVS':<42} {'ClinVar':>7} {'t2g2':>6} {'LV2':>6} {'+ours':>6} {'+lv2':>5}")
    print("-"*80)
    for r in rows:
        flag = " +" if r["only_t2g2"] > 0 else ("  =" if r["overlap"] > 0 else "")
        print(f"{r['gene']:<10} {r['hgvs']:<42} {r['n_clinvar']:>7} "
              f"{r['n_t2g2']:>6} {r['n_lv2']:>6} {r['only_t2g2']:>6} {r['only_lv2']:>5}{flag}")
    print("="*80)

    n = len(rows)
    print(f"\nVariants: {n}")
    for gene in sorted(HHT_GENES):
        sub = [r for r in rows if r["gene"] == gene]
        if sub:
            print(f"  {gene}: {len(sub)} variants")
    print(f"\nPMID totals:")
    print(f"  ClinVar citations:   {sum(r['n_clinvar'] for r in rows)}")
    print(f"  text2gene2 total:    {sum(r['n_t2g2'] for r in rows)}")
    print(f"  LitVar2 total:       {sum(r['n_lv2'] for r in rows)}")
    print(f"\nPer-variant comparison:")
    print(f"  t2g2 > LitVar2:      {sum(1 for r in rows if r['only_t2g2'] > 0)}/{n}")
    print(f"  LitVar2 > t2g2:      {sum(1 for r in rows if r['only_lv2']  > 0)}/{n}")
    print(f"  Equal (both > 0):    {sum(1 for r in rows if r['overlap'] > 0 and r['only_t2g2'] == 0 and r['only_lv2'] == 0)}/{n}")
    print(f"  Both zero:           {sum(1 for r in rows if r['n_t2g2'] == 0 and r['n_lv2'] == 0)}/{n}")


async def main_async(args):
    citation_map  = fetch_citation_counts()
    variants      = fetch_hht_variants(citation_map)
    lv2_lookup    = build_litvar2_lookup(variants)
    t2g2_lookup   = await query_all_t2g2(variants, args.api, args.concurrency)

    rows = []
    for v in variants:
        t2g2_pmids   = t2g2_lookup.get(v["hgvs"], set())
        lv2_pmids    = lv2_lookup.get(v["hgvs"], set())
        clinvar_pmids = v["clinvar_pmids"]
        rows.append({
            **v,
            "n_clinvar":  len(clinvar_pmids),
            "n_t2g2":     len(t2g2_pmids),
            "n_lv2":      len(lv2_pmids),
            "overlap":    len(t2g2_pmids & lv2_pmids),
            "only_t2g2":  len(t2g2_pmids - lv2_pmids),
            "only_lv2":   len(lv2_pmids  - t2g2_pmids),
            "clinvar_pmids": None,   # drop the set before writing
        })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, out)
    print_summary(rows)
    print(f"\nFull results → {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--api", default="http://loki.local",
                    help="text2gene2 base URL (default: http://loki.local)")
    ap.add_argument("--concurrency", type=int, default=8,
                    help="Concurrent text2gene2 requests (default: 8)")
    ap.add_argument("--output", default="benchmarks/hht_benchmark.tsv")
    args = ap.parse_args()
    asyncio.run(main_async(args))


if __name__ == "__main__":
    main()
