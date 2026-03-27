# text2gene2 — Claude Code Project Instructions

## What This Is

text2gene2 is a bioinformatics tool: given an HGVS-formatted genetic variant,
find all relevant PubMed literature by querying multiple sources in parallel.
It is deployed as a web service at text2gene.org (running on loki.local).

## Repository Layout

```
src/text2gene2/
  lvg/          VariantValidator API client — HGVS expansion across transcript versions
  sources/      One async source per external database (litvar2, clinvar, google_cse, europepmc)
  pipeline/     fanout.py (asyncio.gather), validate.py (abstract-level validation)
  api/          FastAPI app + Jinja2 web UI
  cli.py        hgvs2pmid CLI (typer)
  cache.py      Redis + in-memory fallback
  rate_limit.py Token-bucket limiters per external service
benchmarks/     ClinVar-derived variant datasets and comparison scripts
deploy/         nginx config, systemd service, install.sh, update-dns.sh
```

## Infrastructure

### Production
- **App**: loki.local (192.168.0.3), systemd service `text2gene2`, 4 uvicorn workers
- **Redis**: loki.local, default port 6379
- **nginx**: loki.local port 80 → uvicorn 127.0.0.1:8000
- **DNS**: text2gene.org via Gandi LiveDNS API (key in nthmost-systems/.secrets/gandi-apikey)
- **Dynamic DNS cron**: `/usr/local/bin/text2gene-update-dns` every 15 min on loki
- **Deployment**: `rsync -az src/ loki.local:/opt/text2gene2/src/ && ssh loki.local "sudo systemctl restart text2gene2"`

### LLM Router
**spartacus.local:4000** — OpenAI-compatible LiteLLM router across local and cloud models.
Use for any LLM-assisted features (validation, classification, query expansion).

Available models:
| ID | Notes |
|---|---|
| `heavy/reasoning` | → Claude 3.5 Sonnet (cloud) — use for reasoning tasks |
| `heavy/coding` | → Claude 3.5 Sonnet (cloud) — use for code generation |
| `gpt4` | → GPT-4 class (cloud) |
| `loki/qwen3-coder-30b` | Local on loki — strong coder, good for structured output |
| `loki/deepseek-r1` | Local on loki — reasoning model |
| `loki/qwen-coder-14b` | Local on loki |
| `loki/llama3-8b` | Local on loki — fast, small tasks |
| `spartacus/qwen3-14b` | Local on spartacus |
| `spartacus/qwen-coder-14b` | Local on spartacus |
| `local/llama-70b` | Router-picked local 70B |
| `local/mistral-24b` | Router-picked local 24B |

Planned uses in text2gene2:
- **Abstract validation** (`pipeline/validate.py`): replace regex with `heavy/reasoning` prompt
- **Evidence classification**: classify papers as functional/clinical/GWAS/case report/review
- **Query expansion**: generate alternate search terms for zero-result rare variants

Base URL for OpenAI-compatible calls: `http://spartacus.local:4000/v1`
No auth key required on LAN.

## Data Sources

| Source | Purpose | Rate limit |
|---|---|---|
| VariantValidator REST API | LVG expansion across transcript versions | 1.5 req/sec (public) |
| LitVar2 API (NCBI) | Text-mined variant→PMID (replaces PubTator) | ~5 req/sec polite |
| ClinVar via Entrez eutils | Curated variant→PMID associations | 3/sec (10 with NCBI key) |
| Google CSE | Papers outside PubMed entirely | 100/day free |
| Europe PMC | Open-access + preprints | ~5 req/sec |

LitVar2 bulk FTP: `https://ftp.ncbi.nlm.nih.gov/pub/lu/LitVar/litvar2_variants.json.gz`
— 1.8GB, newline-delimited JSON, `gene` field is a list, `all_hgvs` format is `"count|short_hgvs"`

## Known Issues / Watch Out For

- **LitVar2 query guard**: protein notations starting with `*` (stop-codon extensions)
  mangle to `"GENE p."` which matches any top variant for that gene. Fixed via
  `_MEANINGFUL_RE` in `sources/litvar2.py` — don't remove this guard.

- **VariantValidator rate limit**: 2 req/sec on public API. Running the benchmark
  scripts concurrently will trigger 429s. Use `--concurrency 4` or less for bulk runs.

- **Redis cache keys**: LitVar2 cache keyed on `litvar2:{input_hgvs}`. If a bug
  produces wrong results, invalidate with `redis-cli DEL "litvar2:{hgvs}"` on loki.

- **Google CSE**: both engine IDs (`GOOGLE_CSE_WHITELIST_ID`, `GOOGLE_CSE_SCHEMA_ID`)
  are live but need a valid `GOOGLE_API_KEY` in `/opt/text2gene2/.env` on loki.
  The engines were created under naomi@nthmost.com's Google account.

## Benchmarks

- `benchmarks/rare_variants.tsv` — 150 rare-disease P/LP variants (≤3 ClinVar citations)
  across 150 genes. Built from ClinVar + var_citations.txt. The general "low bar" dataset.

- `benchmarks/hht_benchmark.tsv` — 1,586 P/LP HHT variants (ENG/ACVRL1/SMAD4/GDF2).
  Comparison: t2g2=26,465 PMIDs vs LitVar2=2,283 PMIDs. Note: LitVar2 numbers were
  derived from bulk FTP streaming with HGVS short-form matching (no rsID matches —
  most rare HHT variants lack rsIDs). Re-run after any pipeline changes.

- Rebuild either dataset: `python benchmarks/build_*.py` (streams from NCBI FTP, ~5 min)

## Domain / DNS

text2gene.org registered via Gandi (API key in nthmost-systems/.secrets/gandi-apikey).
HTTPS pending — run `~/get-cert.sh` on loki once DNS propagates (check: `dig text2gene.org @8.8.8.8`).
Auto-renew is enabled on Gandi.
