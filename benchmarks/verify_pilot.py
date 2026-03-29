#!/usr/bin/env python3
"""
Pilot verification: fetch literature for 5 rare-disease variants from text2gene2
and LitVar2, then enrich each PMID with metadata + fulltext availability via
metapub. Output a TSV for manual relevance annotation.

Usage:
  pip install -e ".[benchmark]"
  python benchmarks/verify_pilot.py [--api http://loki.local] [--output benchmarks/pilot_report.tsv]

Environment (.env or shell):
  TEXT2GENE2_API   — text2gene2 base URL (default: http://loki.local)
  NCBI_API_KEY     — optional but recommended to avoid rate limiting
"""

import argparse
import asyncio
import csv
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import httpx

from text2gene2.config import settings

try:
    from metapub import PubMedFetcher
except ImportError:
    print("ERROR: metapub not installed. Run: pip install -e '.[benchmark]'", file=sys.stderr)
    sys.exit(1)

# ── Pilot variant list ────────────────────────────────────────────────────────

PILOT_VARIANTS = [
    {"hgvs": "NM_000492.4:c.3963+1G>A", "gene": "CFTR",  "disease": "Cystic fibrosis"},
    {"hgvs": "NM_000070.3:c.257C>T",    "gene": "CAPN3", "disease": "Autosomal recessive limb-girdle muscular dystrophy type 2A"},
    {"hgvs": "NM_000521.4:c.771+5G>C",  "gene": "HEXB",  "disease": "Sandhoff disease (related)"},
    {"hgvs": "NM_020247.5:c.815G>A",    "gene": "COQ8A", "disease": "Autosomal recessive ataxia due to ubiquinone deficiency"},
    {"hgvs": "NM_001876.4:c.1216C>T",   "gene": "CPT1A", "disease": "Carnitine palmitoyl transferase 1A deficiency"},
]

LITVAR2_AUTOCOMPLETE = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/"
LITVAR2_PUBLICATIONS = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/publications/{litvar_id}/"

# ── Data model ────────────────────────────────────────────────────────────────

@dataclass
class ArticleRow:
    hgvs: str
    gene: str
    disease: str
    pmid: int
    source: str          # "text2gene2", "litvar2", or "both"
    title: str = ""
    abstract: str = ""   # first 300 chars
    pmc_id: str = ""
    fulltext_url: str = ""
    relevance: str = ""  # blank — for manual annotation

# ── LitVar2 queries ───────────────────────────────────────────────────────────

async def litvar2_pmids(hgvs: str, client: httpx.AsyncClient) -> set[int]:
    """
    Two-step LitVar2 lookup:
      1. Autocomplete to resolve the HGVS → LitVar2 variant ID
      2. Fetch publications for that ID
    """
    try:
        r = await client.get(LITVAR2_AUTOCOMPLETE, params={"query": hgvs}, timeout=30)
        hits = r.json()
        if not hits:
            return set()
        litvar_id = hits[0].get("_id") or hits[0].get("id")
        if not litvar_id:
            return set()

        r2 = await client.get(
            LITVAR2_PUBLICATIONS.format(litvar_id=litvar_id),
            timeout=30,
        )
        data = r2.json()
        pmids = {int(p) for p in (data.get("pmids") or [])}
        return pmids

    except Exception as e:
        print(f"  LitVar2 error for {hgvs}: {e}", file=sys.stderr)
        return set()


# ── text2gene2 queries ────────────────────────────────────────────────────────

async def t2g2_pmids(hgvs: str, api_base: str, client: httpx.AsyncClient) -> set[int]:
    try:
        r = await client.get(f"{api_base}/api/v2/hgvs2pmid/{hgvs}", timeout=90)
        return {c["pmid"] for c in r.json().get("citations", [])}
    except Exception as e:
        print(f"  text2gene2 error for {hgvs}: {e}", file=sys.stderr)
        return set()


# ── metapub enrichment ────────────────────────────────────────────────────────

def enrich_pmids(pmid_set: set[int], fetcher: PubMedFetcher, delay: float = 0.34) -> dict[int, dict]:
    """
    Fetch article metadata for each PMID via metapub.
    Returns {pmid: {title, abstract, pmc_id, fulltext_url}}.
    Respects NCBI rate limit (~3 req/s without API key, ~10 with).
    """
    results = {}
    total = len(pmid_set)
    for i, pmid in enumerate(sorted(pmid_set), 1):
        print(f"  metapub {i}/{total}: PMID {pmid}", end="\r", flush=True)
        try:
            art = fetcher.article_by_pmid(str(pmid))
            results[pmid] = {
                "title":        art.title or "",
                "abstract":     (art.abstract or "")[:300],
                "pmc_id":       art.pmc or "",
                "fulltext_url": art.url or "",
            }
        except Exception as e:
            print(f"\n  metapub error for PMID {pmid}: {e}", file=sys.stderr)
            results[pmid] = {"title": "", "abstract": "", "pmc_id": "", "fulltext_url": ""}
        time.sleep(delay)
    print()  # clear \r line
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

async def gather_pmids(variants: list[dict], api_base: str) -> dict[str, dict[str, set[int]]]:
    """
    Returns {hgvs: {"t2g2": set, "lv2": set}} for all pilot variants.
    """
    print("Querying text2gene2 and LitVar2 for all pilot variants...", flush=True)
    results: dict[str, dict[str, set[int]]] = {}

    async with httpx.AsyncClient() as client:
        for v in variants:
            hgvs = v["hgvs"]
            print(f"  {hgvs}", flush=True)
            t2g2, lv2 = await asyncio.gather(
                t2g2_pmids(hgvs, api_base, client),
                litvar2_pmids(hgvs, client),
            )
            results[hgvs] = {"t2g2": t2g2, "lv2": lv2}
            print(f"    text2gene2: {len(t2g2)} PMIDs  |  LitVar2: {len(lv2)} PMIDs", flush=True)

    return results


def build_rows(variants: list[dict],
               pmid_map: dict[str, dict[str, set[int]]],
               article_map: dict[int, dict]) -> list[ArticleRow]:
    rows = []
    for v in variants:
        hgvs = v["hgvs"]
        t2g2 = pmid_map[hgvs]["t2g2"]
        lv2  = pmid_map[hgvs]["lv2"]
        all_pmids = t2g2 | lv2

        for pmid in sorted(all_pmids):
            in_t2g2 = pmid in t2g2
            in_lv2  = pmid in lv2
            if in_t2g2 and in_lv2:
                source = "both"
            elif in_t2g2:
                source = "text2gene2"
            else:
                source = "litvar2"

            meta = article_map.get(pmid, {})
            rows.append(ArticleRow(
                hgvs=hgvs,
                gene=v["gene"],
                disease=v["disease"],
                pmid=pmid,
                source=source,
                title=meta.get("title", ""),
                abstract=meta.get("abstract", ""),
                pmc_id=meta.get("pmc_id", ""),
                fulltext_url=meta.get("fulltext_url", ""),
            ))
    return rows


def write_tsv(rows: list[ArticleRow], path: Path) -> None:
    fields = ["hgvs", "gene", "disease", "pmid", "source",
              "title", "abstract", "pmc_id", "fulltext_url", "relevance"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))
    print(f"\nWrote {len(rows)} rows to {path}")


def print_summary(variants: list[dict],
                  pmid_map: dict[str, dict[str, set[int]]]) -> None:
    print("\n" + "="*72)
    print(f"{'GENE':<8} {'t2g2':>6} {'lv2':>6} {'both':>6} {'only t2g2':>10} {'only lv2':>9}")
    print("-"*72)
    for v in variants:
        hgvs = v["hgvs"]
        t2g2 = pmid_map[hgvs]["t2g2"]
        lv2  = pmid_map[hgvs]["lv2"]
        print(f"{v['gene']:<8} {len(t2g2):>6} {len(lv2):>6} "
              f"{len(t2g2 & lv2):>6} {len(t2g2 - lv2):>10} {len(lv2 - t2g2):>9}")
    print("="*72)


async def main_async(args: argparse.Namespace) -> None:
    # Step 1: PMID retrieval
    pmid_map = await gather_pmids(PILOT_VARIANTS, args.api)

    all_pmids: set[int] = set()
    for v in PILOT_VARIANTS:
        all_pmids |= pmid_map[v["hgvs"]]["t2g2"] | pmid_map[v["hgvs"]]["lv2"]

    print(f"\n{len(all_pmids)} unique PMIDs to enrich via metapub...", flush=True)

    # Step 2: metapub enrichment
    # Use faster rate if NCBI_API_KEY is configured
    delay = 0.11 if settings.ncbi_api_key else 0.34
    fetcher = PubMedFetcher()
    article_map = enrich_pmids(all_pmids, fetcher, delay=delay)

    # Step 3: assemble and write
    rows = build_rows(PILOT_VARIANTS, pmid_map, article_map)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, out)
    print_summary(PILOT_VARIANTS, pmid_map)

    n_fulltext = sum(1 for r in rows if r.pmc_id)
    print(f"\nFulltext available (PMC): {n_fulltext}/{len(rows)} articles")
    print(f"\nNext step: open {out} and fill in the 'relevance' column")
    print("  Values: relevant | partial | not_relevant | unsure")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--api", default=settings.text2gene2_api,
                    help=f"text2gene2 base URL (default: {settings.text2gene2_api})")
    ap.add_argument("--output", default="benchmarks/pilot_report.tsv",
                    help="Output TSV path (default: benchmarks/pilot_report.tsv)")
    args = ap.parse_args()
    asyncio.run(main_async(args))


if __name__ == "__main__":
    main()
