# Data Sources & Access Patterns

## Overview

text2gene2 queries multiple data sources to find literature for a variant. Each
source has its own access pattern, rate limits, caching strategy, and failure
mode. This document describes how every external system is accessed, what data
it provides, and how results are cached.

## Architecture

```
User query (HGVS string)
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  LVG Expansion (VariantValidator REST API)               │
│  → LVGResult: gene_symbol, hgvs_c[], hgvs_p[],         │
│    hgvs_g[], rsids[]                                     │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Query Expansion (pipeline/expand.py)                    │
│  → ExpandedVariant: slang forms, 1-letter protein,      │
│    gene synonyms (from HGNC via Postgres)                │
└───────────────────────┬─────────────────────────────────┘
                        │
          ┌─────────────┼──────────────┐──────────────┐
          ▼             ▼              ▼              ▼
     ┌─────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐
     │ LitVar2 │  │ ClinVar  │  │ EuropePMC│  │ Google   │
     │  (API)  │  │  (API)   │  │  (API)   │  │  CSE     │
     └────┬────┘  └────┬─────┘  └────┬─────┘  └────┬─────┘
          └─────────────┴──────────────┴──────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Merge & Deduplicate → CitationTable                     │
│  → Enrich (pubmed.article via Postgres, then NCBI)      │
│  → Validate (abstracts via NCBI efetch)                  │
└─────────────────────────────────────────────────────────┘
```

## External APIs (HTTP)

### VariantValidator

- **Purpose:** LVG expansion — HGVS input → all equivalent representations
  across transcript versions and sequence types
- **URL:** `https://rest.variantvalidator.org/VariantValidator/variantvalidator/{build}/{hgvs}/all`
- **Auth:** None (public API)
- **Rate limit:** 1.5 req/sec (self-imposed; public API asks <= 2/sec)
- **Timeout:** 30s
- **Client:** `httpx.AsyncClient` (new client per request)
- **Cache key:** `lvg:vv:{hgvs}` in Redis, TTL 7 days
- **Failure mode:** Returns minimal `LVGResult` with warning; pipeline continues
  with just the input HGVS
- **Code:** `lvg/expander.py`

### LitVar2

- **Purpose:** Text-mined variant→PMID associations (NCBI)
- **URL:** `https://www.ncbi.nlm.nih.gov/research/litvar2-api`
- **Endpoints used:**
  - `/variant/autocomplete/?query={hgvs}` → resolves HGVS to rsIDs
  - `/variant/get/{varid}/publications` → rsID → list of PMIDs
- **Auth:** None
- **Rate limit:** 5 req/sec (self-imposed; no documented limit)
- **Timeout:** 15s per request
- **Client:** `httpx.AsyncClient` (one client per query, shared across sub-requests)
- **Cache key:** `litvar2:{input_hgvs}` in Redis, TTL 24 hours
- **Failure mode:** Returns empty PMID list; individual autocomplete/publication
  failures are logged and skipped
- **Query strategy:** Tries autocomplete with input HGVS, then gene + protein
  short form, then gene + coding short form. Merges rsIDs from all routes,
  then fetches publications for each rsID.
- **Known issue:** Protein notations starting with `*` (stop-codon extensions)
  mangle to degenerate queries. Guarded by `_MEANINGFUL_RE`.
- **Code:** `sources/litvar2.py`

### ClinVar (via NCBI Entrez eutils)

- **Purpose:** Curated variant→PMID associations
- **URL:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`
- **Endpoints used:**
  - `esearch.fcgi?db=clinvar&term={hgvs}` → ClinVar VariationIDs
  - `elink.fcgi?dbfrom=clinvar&db=pubmed&id={ids}` → PMIDs
- **Auth:** Optional NCBI API key (`NCBI_API_KEY` env var). With key: 10 req/sec.
  Without: 3 req/sec.
- **Rate limit:** 3 req/sec (self-imposed; matches unauthenticated limit)
- **Timeout:** 15s
- **Client:** `httpx.AsyncClient` (one client per query)
- **Cache key:** `clinvar:{input_hgvs}` in Redis, TTL 24 hours
- **Failure mode:** Returns empty PMID list
- **Query strategy:** Builds ClinVar search terms from up to 3 coding HGVS forms,
  2 genomic forms, and 3 rsIDs. Falls back to raw input HGVS if LVG is empty.
- **Code:** `sources/clinvar.py`

### Europe PMC

- **Purpose:** Open-access literature including preprints and gray literature
- **URL:** `https://www.ebi.ac.uk/europepmc/webservices/rest/search`
- **Auth:** None
- **Rate limit:** 5 req/sec (self-imposed; API is generous)
- **Timeout:** 20s
- **Client:** `httpx.AsyncClient` (new client per request)
- **Cache key:** `europepmc:{input_hgvs}` in Redis, TTL 24 hours
- **Failure mode:** Returns empty PMID list with error string
- **Query strategy:** Uses expanded query from `pipeline/expand.py`:
  `("gene" OR "alias1" OR ...) AND ("p.form" OR "1letter" OR "c.form" OR "slang" OR "rsID")`
  See `QUERY_EXPANSION.md` for full details.
- **Response handling:** Requests `resultType=core`, `pageSize=50`, sorted by
  citation count descending. Only results with a numeric `pmid` field are kept.
- **Code:** `sources/europepmc.py`

### Google CSE

- **Purpose:** Papers outside PubMed (publisher sites, institutional repos)
- **URL:** `https://www.googleapis.com/customsearch/v1`
- **Auth:** `GOOGLE_API_KEY` env var (required)
- **Rate limit:** 1 req/sec (burst ok; free tier is 100 queries/day)
- **Timeout:** 20s
- **Client:** `httpx.AsyncClient` (one client per query)
- **Cache key:** `google:{input_hgvs}` in Redis, TTL 7 days
- **Failure mode:** Returns empty PMID list. If no API key configured, returns
  immediately with error.
- **Engines:** Two CSE engine IDs (whitelist + schema). Both are queried.
- **PMID extraction:** Regex on result URLs and snippets for PubMed links.
- **Note:** Not yet wired to query expansion (issue #17). Currently uses simple
  gene + protein/coding form queries.
- **Code:** `sources/google_cse.py`

### NCBI efetch (enrichment + validation)

- **Purpose:** Fetch article metadata (title, authors, journal, year, DOI, abstract)
  for display and abstract-level validation
- **URL:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`
- **Auth:** Optional NCBI API key
- **Rate limit:** Shares the `clinvar` rate limiter (3 req/sec)
- **Timeout:** 30s
- **Batch size:** 100 PMIDs per efetch call (enrichment), 200 per call (validation)
- **Used by:**
  - `pipeline/enrich.py` — fetches metadata for display (title, authors, etc.)
  - `pipeline/validate.py` — fetches abstracts for variant-mention scoring
- **Enrichment is a last resort:** enrich.py tries Redis cache first, then local
  `pubmed.article` DB, then NCBI efetch only for remaining misses.
- **Code:** `pipeline/enrich.py`, `pipeline/validate.py`

## PostgreSQL (medgen on loki.local)

### Connection details

- **Host:** loki.local (192.168.0.3)
- **Port:** 5432
- **Database:** medgen
- **User:** medgen
- **Password:** medgen
- **Auth:** scram-sha-256, open to LAN (192.168.0.0/16 + IPv6)

### Connection management

All Postgres access goes through `db.py`, which provides a single shared
connection via `get_medgen_conn()`:

- **Config:** `settings.medgen_db_url` (env var `MEDGEN_DB_URL`, default
  `postgresql://medgen:medgen@loki.local:5432/medgen`)
- **Connection:** Lazy singleton `psycopg2` with `readonly=True, autocommit=True`
- **Reconnect:** If a query fails, call `reset_medgen_conn()` to force reconnect
  on next access
- **Async:** Synchronous psycopg2 wrapped in `asyncio.to_thread()` by callers.
  Adequate for current query patterns (fast single queries) but would need
  connection pooling for higher concurrency.

All modules that need Postgres import from `db.py`:
- `pipeline/enrich.py` — queries `pubmed.article` for metadata
- `pipeline/expand.py` — queries `hgnc.gene` for gene synonyms
- Future: LOVD, local ClinVar, PubTator sources

### Schemas and tables currently accessed

| Schema.Table | Used by | Purpose | Rows |
|-------------|---------|---------|------|
| `hgnc.gene` | `pipeline/expand.py` | Gene synonym lookup (alias_symbols, prev_symbols) | ~43K |
| `pubmed.article` | `pipeline/enrich.py` | Article metadata (title, authors, journal, year, DOI, abstract) | ~40M |

### Schemas and tables available but not yet used

| Schema.Table | Potential use | Rows | Issue |
|-------------|--------------|------|-------|
| `lovd.lsdb` | LOVD instance lookup (gene → LOVD URLs) | 84,660 | #6 |
| `clinvar.var_citations` | Local ClinVar PMID lookup (bypass eutils API) | 4.5M | #8 |
| `clinvar.variant_summary` | HGVS → VariationID resolution for local ClinVar | 374K P/LP | #8 |
| `gene.pubmed` | Gene-level PMID confirmation signal | 2.2M | #9 |
| `medgen.pubmed` | Disease concept → PMID discovery | 47.8M | #10 |
| `pubtator.mutation_mention` | Text-mined mutation→PMID (local PubTator) | 7.7M | #7, #12 |
| `pubtator.gene_mention` | Gene mention → PMID | 83.2M | — |
| `hpo.gene_phenotype` | Gene → HPO phenotype terms | — | #10 |
| `orphanet.gene_association` | Gene → rare disease associations | — | #10 |
| `pmc.article` | Full-text article content for validation | 5.1M | #11 |

## Redis (loki.local)

- **URL:** `redis://localhost:6379/0` (from app's perspective on loki)
- **Client:** `redis.asyncio` (async Redis client)
- **Fallback:** In-memory dict if Redis is unavailable (dev mode)
- **Connection:** Created and pinged on every cache operation. No connection
  pooling. This works but is not ideal for high-throughput use.

### Cache key scheme

| Prefix | Source | TTL | Example |
|--------|--------|-----|---------|
| `lvg:vv:{hgvs}` | VariantValidator LVG result | 7 days | `lvg:vv:NM_000020.3:c.1120C>T` |
| `litvar2:{hgvs}` | LitVar2 PMID list | 24 hours | `litvar2:NM_000020.3:c.1120C>T` |
| `clinvar:{hgvs}` | ClinVar PMID list | 24 hours | `clinvar:NM_000020.3:c.1120C>T` |
| `europepmc:{hgvs}` | EuropePMC PMID list | 24 hours | `europepmc:NM_000020.3:c.1120C>T` |
| `google:{hgvs}` | Google CSE PMID list | 7 days | `google:NM_000020.3:c.1120C>T` |
| `meta:{pmid}` | Article metadata (title, authors, etc.) | 30 days | `meta:12345678` |

All values are JSON-serialized. Cache keys use the **input HGVS** (not expanded
forms), so a cache hit short-circuits the entire source query.

**Cache invalidation:** Manual only. If a bug produces wrong cached results:
```bash
redis-cli -h loki.local DEL "litvar2:NM_000020.3:c.1120C>T"
```

### In-memory cache

`pipeline/expand.py` maintains its own in-memory cache for gene synonyms
(`_synonym_cache` dict). This never expires within a process lifetime, which is
fine because gene names change rarely and workers restart on deploy.

## Rate Limiting

All external APIs are rate-limited via `rate_limit.py` using a simple async
token-bucket implementation. One `RateLimiter` instance per service, shared
across all requests within a uvicorn worker process.

| Service | Rate | Notes |
|---------|------|-------|
| VariantValidator | 1.5 req/sec | Public API requests <= 2/sec |
| LitVar2 | 5 req/sec | No documented limit; polite default |
| ClinVar / NCBI eutils | 3 req/sec | 10/sec with API key (not yet configured in prod) |
| Google CSE | 1 req/sec | 100/day free tier; can burst within that |
| Europe PMC | 5 req/sec | Generous public API |

Rate limiters are per-worker. With 4 uvicorn workers, the effective rate to each
service can be up to 4x the configured rate. For well-behaved external services
this is fine; for VariantValidator (strict 2/sec) this could be a problem under
heavy concurrent load.

## Pipeline Flow

For a single variant query, the data access sequence is:

1. **Redis cache check** for LVG (`lvg:vv:{hgvs}`)
2. If miss: **VariantValidator API** call → cache result
3. **Query expansion** using LVG result:
   - **Postgres** (`hgnc.gene`) for gene synonyms (or in-memory cache)
   - Deterministic string transforms for slang forms (no DB/API)
4. **Parallel fan-out** to all 4 sources (each checks Redis first):
   - LitVar2: Redis → autocomplete API → publications API → cache
   - ClinVar: Redis → esearch API → elink API → cache
   - EuropePMC: Redis → search API → cache
   - Google CSE: Redis → search API → cache
5. **Merge & deduplicate** PMIDs into CitationTable
6. **Enrich** metadata (for web UI):
   - Redis (`meta:{pmid}`) → Postgres (`pubmed.article`) → NCBI efetch
7. **Validate** (optional): NCBI efetch for abstracts → regex scoring
