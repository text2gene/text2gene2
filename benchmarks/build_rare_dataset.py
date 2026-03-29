#!/usr/bin/env python3
"""
Build a benchmark dataset of rare-disease HGVS variants with minimal existing literature.

Sources:
  ClinVar variant_summary.txt.gz  — variant metadata, Orphanet disease links
  ClinVar var_citations.txt       — AlleleID → PMID associations

Strategy:
  1. Build citation count map from var_citations.txt (AlleleID → n PMIDs)
  2. Stream variant_summary.txt for GRCh38 P/LP variants linked to Orphanet diseases
  3. Filter for ≤ MAX_CITATIONS citations (the "low bar")
  4. Parse clean HGVS from the Name column
  5. Sample a diverse set capped at TARGET_SIZE variants
  6. Write benchmarks/rare_variants.tsv

Usage:
  python benchmarks/build_rare_dataset.py [--max-citations 3] [--target 150] [--output rare_variants.tsv]

Environment (.env or shell):
  NCBI_API_KEY   — not used directly here, but loaded for consistency with other benchmark scripts
"""

import argparse
import csv
import gzip
import io
import random
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

from text2gene2.config import settings  # noqa: F401  (loads .env)

VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
VAR_CITATIONS_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"

# Regex to pull a clean transcript HGVS out of ClinVar's Name field.
# Name looks like: NM_014855.3(AP5Z1):c.80_83delinsTGCT(p.Ala27LeufsTer123)
# We want:         NM_014855.3:c.80_83delinsTGCT
_HGVS_RE = re.compile(r"(N[MR]_\d+\.\d+)(?:\([^)]+\))?:(c\.[^\s(,]+)")


def parse_hgvs(name: str) -> str | None:
    m = _HGVS_RE.search(name)
    return f"{m.group(1)}:{m.group(2)}" if m else None


def stream_url(url: str):
    """Yield lines from a (possibly gzipped) URL without loading into memory."""
    print(f"  Fetching {url} ...", flush=True)
    with urllib.request.urlopen(url) as resp:
        raw = resp.read()
    if url.endswith(".gz"):
        f = gzip.open(io.BytesIO(raw), "rt", encoding="utf-8", errors="replace")
    else:
        f = io.StringIO(raw.decode("utf-8", errors="replace"))
    yield from f


def build_citation_counts(max_citations: int) -> dict[int, int]:
    """
    Returns {allele_id: pmid_count} for alleles with <= max_citations PMIDs.
    Alleles not in var_citations.txt have 0 citations.
    """
    print("Step 1: counting citations per allele...", flush=True)
    counts: dict[int, int] = defaultdict(int)
    for line in stream_url(VAR_CITATIONS_URL):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        allele_id_str = parts[0]
        citation_source = parts[4] if len(parts) > 4 else ""
        if citation_source == "PubMed":
            try:
                counts[int(allele_id_str)] += 1
            except ValueError:
                pass
    print(f"  {len(counts):,} alleles have ≥1 PubMed citation", flush=True)
    return counts


def collect_variants(citation_counts: dict[int, int], max_citations: int) -> list[dict]:
    """
    Stream variant_summary and collect GRCh38 P/LP Orphanet-linked variants
    with citation count <= max_citations.
    """
    print("Step 2: scanning variant_summary for rare-disease P/LP variants...", flush=True)
    candidates: list[dict] = []
    seen_hgvs: set[str] = set()
    n_rows = 0

    for line in stream_url(VARIANT_SUMMARY_URL):
        if line.startswith("#AlleleID") or line.startswith("AlleleID"):
            # header — extract column indices
            cols = line.lstrip("#").rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(cols)}
            continue

        n_rows += 1
        if n_rows % 200_000 == 0:
            print(f"  ...{n_rows:,} rows scanned, {len(candidates):,} candidates so far", flush=True)

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 20:
            continue

        def get(col):
            i = idx.get(col)
            return parts[i] if i is not None and i < len(parts) else ""

        # Only GRCh38, only Pathogenic/LP
        if get("Assembly") != "GRCh38":
            continue
        if get("ClinSigSimple") != "1":
            continue

        # Must be linked to an Orphanet disease
        phenotype_ids = get("PhenotypeIDS")   # note: IDS not IDs in actual file
        if "Orphanet:" not in phenotype_ids:
            continue

        # Citation filter
        try:
            allele_id = int(get("AlleleID"))
        except ValueError:
            continue
        n_citations = citation_counts.get(allele_id, 0)
        if n_citations > max_citations:
            continue

        # Parse HGVS
        name = get("Name")
        hgvs = parse_hgvs(name)
        if not hgvs or hgvs in seen_hgvs:
            continue
        seen_hgvs.add(hgvs)

        gene = get("GeneSymbol").split(";")[0].split(",")[0].strip()
        disease = get("PhenotypeList").split("|")[0].strip()
        rsid = get("RS# (dbSNP)").strip()

        candidates.append({
            "allele_id":    allele_id,
            "variation_id": get("VariationID"),
            "hgvs":         hgvs,
            "gene":         gene,
            "disease":      disease,
            "rsid":         rsid if rsid != "-1" else "",
            "n_citations":  n_citations,
            "review_status": get("ReviewStatus"),
        })

    print(f"  Scanned {n_rows:,} rows total. {len(candidates):,} candidates pass filters.", flush=True)
    return candidates


def sample_diverse(candidates: list[dict], target: int) -> list[dict]:
    """
    Sample up to `target` variants, maximising gene and disease diversity.
    Prefer variants with 1-3 citations over 0 (at least some signal exists).
    """
    # Sort: prefer 1-3 citations (some signal) over 0 (possibly just never studied)
    # Within that, shuffle for diversity
    random.seed(42)
    with_signal    = [v for v in candidates if 1 <= v["n_citations"] <= 3]
    without_signal = [v for v in candidates if v["n_citations"] == 0]
    random.shuffle(with_signal)
    random.shuffle(without_signal)
    pool = with_signal + without_signal

    # Greedy gene-diversity sampling
    seen_genes: set[str] = set()
    first_pass, second_pass = [], []
    for v in pool:
        if v["gene"] not in seen_genes:
            first_pass.append(v)
            seen_genes.add(v["gene"])
        else:
            second_pass.append(v)

    selected = (first_pass + second_pass)[:target]
    selected.sort(key=lambda v: (v["gene"], v["hgvs"]))
    return selected


def write_tsv(variants: list[dict], path: Path) -> None:
    fields = ["hgvs", "gene", "disease", "n_citations", "rsid", "allele_id", "variation_id", "review_status"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(variants)
    print(f"\nWrote {len(variants)} variants to {path}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-citations", type=int, default=3,
                    help="Maximum existing PubMed citations (default: 3)")
    ap.add_argument("--target", type=int, default=150,
                    help="Target dataset size (default: 150)")
    ap.add_argument("--output", default="benchmarks/rare_variants.tsv",
                    help="Output TSV path (default: benchmarks/rare_variants.tsv)")
    args = ap.parse_args()

    print(f"Building rare-disease variant benchmark (max_citations={args.max_citations}, target={args.target})")

    citation_counts = build_citation_counts(args.max_citations)
    candidates = collect_variants(citation_counts, args.max_citations)

    if not candidates:
        print("ERROR: no candidates found. Check network access to NCBI FTP.", file=sys.stderr)
        sys.exit(1)

    selected = sample_diverse(candidates, args.target)
    print(f"\nSelected {len(selected)} variants across {len(set(v['gene'] for v in selected))} genes")

    # Summary
    from collections import Counter
    cit_dist = Counter(v["n_citations"] for v in selected)
    for k in sorted(cit_dist):
        print(f"  {k} citations: {cit_dist[k]} variants")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(selected, out)


if __name__ == "__main__":
    main()
