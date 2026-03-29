#!/usr/bin/env python3
"""
Fetch PMC fulltext for articles in pilot_report.tsv.

Uses NCBI efetch to download full-text XML for each article with a PMC ID,
strips markup to plain text, and saves to benchmarks/fulltext/{pmid}.txt.

Articles without a PMC ID are skipped (abstract only available).

Usage:
  python benchmarks/fetch_fulltext.py [--report benchmarks/pilot_report.tsv] [--out benchmarks/fulltext]

Environment (.env or shell):
  NCBI_API_KEY   — optional but recommended (10 req/sec vs 3 req/sec)
"""

import argparse
import csv
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from xml.etree import ElementTree as ET

from text2gene2.config import settings

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_pmc_xml(pmc_id: str, ncbi_key: str, retries: int = 4) -> str | None:
    params = f"db=pmc&id={pmc_id}&rettype=full&retmode=xml"
    if ncbi_key:
        params += f"&api_key={ncbi_key}"
    url = f"{EFETCH_URL}?{params}"
    wait = 2.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None  # not in PMC OA
            if e.code == 429 or e.code >= 500:
                if attempt < retries - 1:
                    time.sleep(wait)
                    wait *= 2
                    continue
            raise
        except OSError:
            # RemoteDisconnected, timeout, etc.
            if attempt < retries - 1:
                time.sleep(wait)
                wait *= 2
                continue
            raise
    return None


def xml_to_text(xml: str) -> str:
    """
    Extract readable body text from PMC full-text XML.
    Pulls article title, abstract, and body paragraphs; strips all tags.
    """
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        # Fall back to naive tag stripping
        return re.sub(r"<[^>]+>", " ", xml)

    chunks = []

    # Title
    for el in root.iter("article-title"):
        if el.text:
            chunks.append(el.text.strip())
        break

    # Abstract
    for el in root.iter("abstract"):
        chunks.append(_el_text(el))
        break

    # Body sections
    for el in root.iter("body"):
        chunks.append(_el_text(el))
        break

    return "\n\n".join(filter(None, chunks))


def _el_text(el) -> str:
    """Concatenate all text content under an element."""
    return " ".join(
        (t.strip() for t in el.itertext() if t.strip()),

    )


def load_report(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--report", default="benchmarks/pilot_report.tsv")
    ap.add_argument("--out", default="benchmarks/fulltext",
                    help="Directory to write fulltext files (default: benchmarks/fulltext)")
    args = ap.parse_args()

    report_path = Path(args.report)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_report(report_path)
    ncbi_key = settings.ncbi_api_key
    delay = 0.11 if ncbi_key else 0.34

    # Deduplicate: one fetch per PMC ID
    seen: set[str] = set()
    to_fetch = []
    for row in rows:
        pmc_id = row.get("pmc_id", "").strip()
        if pmc_id and pmc_id not in seen:
            seen.add(pmc_id)
            to_fetch.append((row["pmid"], pmc_id))

    no_pmc = sum(1 for r in rows if not r.get("pmc_id", "").strip())
    print(f"{len(to_fetch)} articles with PMC IDs to fetch")
    print(f"{no_pmc} articles have no PMC ID (abstract only — skipped)")
    if ncbi_key:
        print("NCBI_API_KEY set — using 10 req/sec rate")
    else:
        print("No NCBI_API_KEY — using 3 req/sec rate")
    print()

    skipped = ok = failed = 0
    for i, (pmid, pmc_id) in enumerate(to_fetch, 1):
        out_path = out_dir / f"{pmid}.txt"
        if out_path.exists():
            skipped += 1
            print(f"  [{i}/{len(to_fetch)}] PMID {pmid} (PMC{pmc_id}) — already exists, skipping")
            continue

        print(f"  [{i}/{len(to_fetch)}] PMID {pmid} (PMC{pmc_id})", end=" ", flush=True)
        xml = fetch_pmc_xml(pmc_id, ncbi_key)
        if xml is None:
            print("— not in PMC OA")
            failed += 1
        else:
            text = xml_to_text(xml)
            out_path.write_text(text, encoding="utf-8")
            words = len(text.split())
            print(f"— {words:,} words")
            ok += 1

        time.sleep(delay)

    print(f"\nDone. {ok} fetched, {skipped} skipped (cached), {failed} not in PMC OA")
    print(f"Fulltext saved to: {out_dir}/")


if __name__ == "__main__":
    main()
