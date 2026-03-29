#!/usr/bin/env python3
"""
Genetic curation of pilot corpus using 3 LLM curators.

Each article in pilot_report.tsv is evaluated by three independent LLMs acting
as genetic curators. They judge whether the paper contains meaningful findings
about the specific variant, not just the gene or disease in general.

Curators:
  - heavy/reasoning   (Claude Sonnet — cloud)
  - loki/deepseek-r1  (DeepSeek R1 — local reasoning model)
  - spartacus/qwen3-14b (Qwen3 — local, independent architecture)

Verdicts:
  relevant     — paper specifically discusses this variant with meaningful findings
                 (functional data, clinical observations, segregation, case report)
  partial      — variant mentioned but without specific analysis
                 (appears in a table, passing mention alongside many others)
  catalog      — variant appears in a large screening list with no specific findings
                 (panel results, population cohort, no per-variant conclusions)
  not_relevant — paper does not appear to involve this specific variant

Usage:
  python benchmarks/curate_pilot.py [--report benchmarks/pilot_report.tsv]
                                     [--fulltext benchmarks/fulltext]
                                     [--out benchmarks/curation_report.tsv]
                                     [--concurrency 3]

Environment (.env or shell):
  TEXT2GENE2_API  — not used here
"""

import argparse
import asyncio
import csv
import json
import os
import re
import sys
from pathlib import Path

import httpx

# Load settings first, then export UTA env vars before metavariant imports
# (metavariant connects to UTA at module import time via get_uta_connection())
from text2gene2.config import settings
os.environ.setdefault("UTA_HOST",   settings.uta_host)
os.environ.setdefault("UTA_PORT",   str(settings.uta_port))
os.environ.setdefault("UTA_USER",   settings.uta_user)
os.environ.setdefault("UTA_PASS",   settings.uta_pass)
os.environ.setdefault("UTA_SCHEMA", settings.uta_schema)

from metavariant import VariantLVG  # noqa: E402 — must come after UTA env setup

LLM_BASE = "http://spartacus.local:4000/v1"

CURATORS = [
    "heavy/reasoning",          # Claude Sonnet — cloud, fast
    "loki/qwen3-coder-30b",     # Qwen3 30B — local, strong, no thinking chain
    "spartacus/qwen3-14b",      # Qwen3 14B — local, independent
]

LEGACY_NAME_PROMPT = """\
You are a specialist in genetic variant nomenclature.

Given a variant in standard HGVS notation, list ALL alternative names, legacy \
notations, and synonyms that might appear in published literature — especially \
older papers predating HGVS standardisation (pre-2000s). Include:
  - Pre-HGVS / exon-based / cDNA numbering variants (e.g. old CFTR numbering \
where 3849+1G>A = current c.3963+1G>A)
  - IVS (intervening sequence) notation for splice-site variants
  - Protein-level trivial names (e.g. ΔF508, p.Phe508del)
  - dbSNP rsID if widely cited
  - Any common shorthand used in the literature for THIS specific gene

Return ONLY a JSON array of strings. No explanation, no markdown.
Example: ["3849+1G>A", "IVS19+1G>A", "rs121908769"]
If you are uncertain or know of no widely-used alternatives, return an empty array [].
"""


async def fetch_legacy_names(hgvs: str, gene: str,
                              client: httpx.AsyncClient) -> list[str]:
    """Ask the LLM to generate legacy/alternative names for this variant."""
    try:
        r = await client.post(
            f"{LLM_BASE}/chat/completions",
            json={
                "model": "heavy/reasoning",
                "messages": [
                    {"role": "system", "content": LEGACY_NAME_PROMPT},
                    {"role": "user",   "content": f"Gene: {gene}\nHGVS: {hgvs}"},
                ],
                "temperature": 0.0,
                "max_tokens": 256,
            },
            timeout=30,
        )
        raw = r.json()["choices"][0]["message"]["content"].strip()
        raw = re.sub(r"<think>.*?</think>", "", raw, flags=re.DOTALL).strip()
        raw = re.sub(r"^```(?:json)?\s*", "", raw, flags=re.IGNORECASE)
        raw = re.sub(r"\s*```$", "", raw)
        names = json.loads(raw)
        return [str(n).strip() for n in names if str(n).strip()]
    except Exception:
        return []


def expand_hgvs_synonyms(hgvs: str) -> dict:
    """
    Use metavariant.VariantLVG (backed by UTA) to expand an HGVS string into
    all known transcript-level synonyms.  Falls back gracefully on error.
    UTA server is configured via .env (UTA_HOST, UTA_PORT, UTA_USER, etc.)
    and exported to os.environ before metavariant is imported.
    """
    try:
        lvg = VariantLVG(hgvs)
        return {
            "hgvs_c": lvg.hgvs_c or [],
            "hgvs_g": lvg.hgvs_g or [],
            "hgvs_p": lvg.hgvs_p or [],
            "hgvs_n": lvg.hgvs_n or [],
        }
    except Exception as e:
        print(f"  [metavariant fallback for {hgvs}: {e}]", flush=True)
        return {"hgvs_c": [hgvs], "hgvs_g": [], "hgvs_p": [], "hgvs_n": []}


def build_synonym_lines(hgvs: str, lvg: dict,
                        clinvar_names: list[str] | None = None) -> list[str]:
    """
    Produce a flat list of all known aliases for the variant so curators know
    what strings to look for in paper text.  Includes:
      - Standard HGVS (coding, genomic, protein) from LVG
      - Short forms without transcript prefix
      - rsIDs
      - ClinVar Name / OtherName fields (captures legacy and trivial notation)
    """
    seen: set[str] = set()
    lines: list[str] = []

    def add(label: str, values):
        for v in (values if isinstance(values, list) else [values]):
            v = str(v).strip()
            if v and v not in seen:
                seen.add(v)
                lines.append(f"  {label}: {v}")
                # Also add the short form without transcript (e.g. "c.3963+1G>A")
                if ":" in v:
                    short = v.split(":", 1)[1]
                    if short not in seen:
                        seen.add(short)
                        lines.append(f"  short form: {short}")

    add("HGVS (coding)", lvg.get("hgvs_c") or [hgvs])
    add("HGVS (genomic)", lvg.get("hgvs_g") or [])
    add("HGVS (protein)", lvg.get("hgvs_p") or [])
    add("HGVS (non-coding)", lvg.get("hgvs_n") or [])
    for rsid in (lvg.get("rsids") or []):
        rs = str(rsid).strip()
        if rs and rs not in seen:
            seen.add(rs)
            lines.append(f"  rsID: rs{rs}" if not rs.startswith("rs") else f"  rsID: {rs}")
    for name in (clinvar_names or []):
        if name not in seen:
            seen.add(name)
            lines.append(f"  legacy/common name: {name}")

    return lines


SYSTEM_PROMPT = """\
You are a genetic curator specialising in variant–disease associations. \
Your task is to evaluate whether a scientific paper contains meaningful \
information about a specific genetic variant.

You will receive:
  - TARGET VARIANT: the HGVS string, gene symbol, and associated disease
  - PAPER: title and text (may be truncated)

Evaluate carefully and return a JSON object with exactly these fields:

{
  "verdict": "<relevant|partial|catalog|not_relevant>",
  "confidence": "<high|medium|low>",
  "reasoning": "<1–2 sentences>"
}

Verdict definitions:
  relevant     — The paper specifically discusses THIS variant with meaningful \
findings: functional assays, clinical/phenotypic data, segregation analysis, \
or a case report centred on this mutation.
  partial      — The variant is mentioned but not specifically analysed: it \
appears in a supplementary table, is listed among many variants screened, or \
is cited in passing without dedicated findings.
  catalog      — The paper is a large variant list or panel screen. The variant \
may appear in a table but no per-variant conclusions are drawn (e.g. "we \
identified 142 variants including…").
  not_relevant — The paper concerns the gene, protein, or disease but does not \
appear to involve this specific variant at all.

Rules:
  - A paper about the gene broadly is NOT relevant unless this variant is discussed.
  - The variant may appear under ANY of the listed synonyms — check all of them.
  - Legacy or short-form notation (e.g. "c.3963+1G>A" without a transcript prefix) \
counts as a match.
  - Favour 'catalog' over 'partial' when the variant is one of dozens with no \
specific commentary.
  - Return ONLY the JSON object. No preamble, no explanation outside the JSON.
"""


def build_user_message(hgvs: str, gene: str, disease: str,
                       synonym_lines: list[str], text: str) -> str:
    synonyms_block = "\n".join(synonym_lines) if synonym_lines else f"  {hgvs}"
    # /no_think suppresses Qwen3's internal reasoning chain, reducing latency
    return (
        f"/no_think\n\n"
        f"TARGET VARIANT\n"
        f"  Gene:    {gene}\n"
        f"  Disease: {disease}\n\n"
        f"KNOWN SYNONYMS FOR THIS VARIANT (any of these may appear in the paper):\n"
        f"{synonyms_block}\n\n"
        f"PAPER\n{text}"
    )


def load_fulltext(pmid: str, fulltext_dir: Path, max_words: int = 1200) -> str:
    path = fulltext_dir / f"{pmid}.txt"
    if not path.exists():
        return ""
    words = path.read_text(encoding="utf-8").split()
    truncated = " ".join(words[:max_words])
    if len(words) > max_words:
        truncated += "\n[… truncated …]"
    return truncated


def load_report(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def parse_verdict(raw: str) -> dict:
    """Extract JSON from model response, tolerating markdown fences and <think> blocks."""
    # Strip <think>...</think> reasoning blocks (DeepSeek-R1, Qwen3)
    raw = re.sub(r"<think>.*?</think>", "", raw, flags=re.DOTALL).strip()
    # Strip markdown code fences
    raw = re.sub(r"^```(?:json)?\s*", "", raw.strip(), flags=re.IGNORECASE)
    raw = re.sub(r"\s*```$", "", raw.strip())
    # Find first {...} block
    m = re.search(r"\{.*\}", raw, re.DOTALL)
    if not m:
        return {"verdict": "error", "confidence": "low", "reasoning": raw[:200]}
    try:
        data = json.loads(m.group())
        valid_verdicts = {"relevant", "partial", "catalog", "not_relevant", "error"}
        if data.get("verdict") not in valid_verdicts:
            data["verdict"] = "error"
        return data
    except json.JSONDecodeError:
        return {"verdict": "error", "confidence": "low", "reasoning": raw[:200]}


async def call_curator(
    model: str,
    user_msg: str,
    client: httpx.AsyncClient,
    sem: asyncio.Semaphore,
) -> dict:
    async with sem:
        try:
            r = await client.post(
                f"{LLM_BASE}/chat/completions",
                json={
                    "model": model,
                    "messages": [
                        {"role": "system", "content": SYSTEM_PROMPT},
                        {"role": "user",   "content": user_msg},
                    ],
                    "temperature": 0.1,
                    "max_tokens": 1024,
                },
                timeout=180,
            )
            r.raise_for_status()
            content = r.json()["choices"][0]["message"]["content"]
            return parse_verdict(content)
        except Exception as e:
            return {"verdict": "error", "confidence": "low", "reasoning": str(e)[:200]}


def majority_verdict(verdicts: list[str]) -> str:
    counts: dict[str, int] = {}
    for v in verdicts:
        counts[v] = counts.get(v, 0) + 1
    top = max(counts, key=lambda k: counts[k])
    if counts[top] >= 2:
        return top
    return "disputed"


async def curate_all(
    rows: list[dict],
    fulltext_dir: Path,
    concurrency: int,
) -> list[dict]:
    sem = asyncio.Semaphore(concurrency)
    results = []
    total = len(rows)

    limits = httpx.Limits(max_connections=concurrency * len(CURATORS) + 4,
                          max_keepalive_connections=concurrency * len(CURATORS) + 4)
    async with httpx.AsyncClient(limits=limits) as client:

        # Pre-fetch LVG synonyms for each unique HGVS
        unique_hgvs = list({row["hgvs"] for row in rows})
        print(f"Pre-fetching synonyms for {len(unique_hgvs)} unique variants...",
              flush=True)
        lvg_map: dict[str, list[str]] = {}
        # Determine gene for each unique hgvs (from row data)
        hgvs_to_gene = {row["hgvs"]: row["gene"] for row in rows}
        for hgvs in unique_hgvs:
            gene         = hgvs_to_gene.get(hgvs, "")
            lvg          = expand_hgvs_synonyms(hgvs)
            legacy_names = await fetch_legacy_names(hgvs, gene, client)
            syn_lines    = build_synonym_lines(hgvs, lvg, legacy_names)
            lvg_map[hgvs] = syn_lines
            print(f"  {hgvs} → {len(syn_lines)} entries "
                  f"(legacy: {', '.join(legacy_names[:4]) or 'none'})", flush=True)
        print()

        for i, row in enumerate(rows, 1):
            pmid    = row["pmid"]
            hgvs    = row["hgvs"]
            gene    = row["gene"]
            disease = row["disease"]

            fulltext = load_fulltext(pmid, fulltext_dir)
            if fulltext:
                paper_text = fulltext
            else:
                parts = []
                if row.get("title"):
                    parts.append(f"Title: {row['title']}")
                if row.get("abstract"):
                    parts.append(f"Abstract: {row['abstract']}")
                paper_text = "\n".join(parts) if parts else "(no text available)"

            user_msg = build_user_message(hgvs, gene, disease,
                                          lvg_map.get(hgvs, []), paper_text)

            print(f"  [{i}/{total}] PMID {pmid} ({gene}) — querying 3 curators...",
                  end=" ", flush=True)

            curator_results = await asyncio.gather(*[
                call_curator(model, user_msg, client, sem)
                for model in CURATORS
            ])

            verdicts = [r.get("verdict", "error") for r in curator_results]
            consensus = majority_verdict(verdicts)
            print(f"{' | '.join(verdicts)} → {consensus}")

            result = {
                "hgvs":       hgvs,
                "gene":       gene,
                "disease":    disease,
                "pmid":       pmid,
                "source":     row.get("source", ""),
                "title":      row.get("title", ""),
                "consensus":  consensus,
            }
            for j, (model, cr) in enumerate(zip(CURATORS, curator_results), 1):
                short = model.split("/")[-1]
                result[f"curator{j}_model"]      = model
                result[f"curator{j}_verdict"]    = cr.get("verdict", "error")
                result[f"curator{j}_confidence"] = cr.get("confidence", "")
                result[f"curator{j}_reasoning"]  = cr.get("reasoning", "")

            results.append(result)

    return results


def write_tsv(results: list[dict], path: Path) -> None:
    if not results:
        return
    fields = list(results[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(results)


def print_summary(results: list[dict]) -> None:
    from collections import Counter
    print("\n" + "="*60)
    print("CONSENSUS VERDICTS")
    print("-"*60)
    by_gene: dict[str, list] = {}
    for r in results:
        by_gene.setdefault(r["gene"], []).append(r["consensus"])

    for gene, verdicts in sorted(by_gene.items()):
        counts = Counter(verdicts)
        parts = "  ".join(f"{v}:{n}" for v, n in sorted(counts.items()))
        print(f"  {gene:<8}  {parts}")

    print("-"*60)
    total_counts = Counter(r["consensus"] for r in results)
    for verdict, n in sorted(total_counts.items()):
        print(f"  {verdict:<15} {n:>4}  ({100*n/len(results):.0f}%)")
    print(f"  {'TOTAL':<15} {len(results):>4}")
    print("="*60)


async def main_async(args: argparse.Namespace) -> None:
    report_path   = Path(args.report)
    fulltext_dir  = Path(args.fulltext)
    out_path      = Path(args.out)

    rows = load_report(report_path)
    print(f"Loaded {len(rows)} rows from {report_path}")
    print(f"Curators: {', '.join(CURATORS)}")
    print(f"Concurrency per curator: {args.concurrency}\n")

    results = await curate_all(rows, fulltext_dir, args.concurrency)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(results, out_path)
    print(f"\nWrote {len(results)} rows to {out_path}")
    print_summary(results)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--report",     default="benchmarks/pilot_report.tsv")
    ap.add_argument("--fulltext",   default="benchmarks/fulltext")
    ap.add_argument("--out",        default="benchmarks/curation_report.tsv")
    ap.add_argument("--concurrency", type=int, default=3,
                    help="Concurrent requests per curator model (default: 3)")
    args = ap.parse_args()
    asyncio.run(main_async(args))


if __name__ == "__main__":
    main()
