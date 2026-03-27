"""
Google Programmable Search Engine source.

Queries two engines:
  - whitelist: cx = GOOGLE_CSE_WHITELIST_ID  (curated publisher list)
  - schema:    cx = GOOGLE_CSE_SCHEMA_ID     (ScholarlyArticle schema)

Returns PMIDs extracted from result URLs (PubMed links) or DOIs resolved to PMIDs.

Free tier: 100 queries/day.  Paid: $5/1000.
We are conservative — one combined query per variant, not per HGVS form.

Google API docs: https://developers.google.com/custom-search/v1/reference/rest/v1/cse/list
"""
import logging
import re

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_CSE_URL = "https://www.googleapis.com/customsearch/v1"
_PMID_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
_PMID_PARAM_RE = re.compile(r"[?&](?:term|id)=(\d+)")


def _extract_pmids_from_results(items: list[dict]) -> list[int]:
    """Pull PubMed PMIDs out of Google CSE result URLs and snippets."""
    pmids: set[int] = set()
    for item in items:
        for field in (item.get("link", ""), item.get("formattedUrl", ""), item.get("htmlSnippet", "")):
            for m in _PMID_RE.finditer(field):
                pmids.add(int(m.group(1)))
            for m in _PMID_PARAM_RE.finditer(field):
                pmids.add(int(m.group(1)))
    return sorted(pmids)


def _build_query(lvg: LVGResult) -> str:
    """
    Build a Google search query from the LVG.

    Google has a 32-term limit. We use gene + the most specific variant forms.
    """
    parts = []
    if lvg.gene_symbol:
        parts.append(lvg.gene_symbol)

    # Prefer protein notation (most human-readable in literature)
    for p in lvg.hgvs_p[:2]:
        parts.append(f'"{p}"')

    # Add coding form if no protein (e.g. intronic)
    if not lvg.hgvs_p:
        for c in lvg.hgvs_c[:2]:
            parts.append(f'"{c}"')

    # rsID if we have one
    for rsid in lvg.rsids[:1]:
        parts.append(rsid)

    # Fall back to raw input
    if not parts:
        parts.append(f'"{lvg.input_hgvs}"')

    return " ".join(parts[:31])  # stay under 32-term limit


async def _cse_query(query: str, cx: str, client: httpx.AsyncClient) -> list[dict]:
    """Run a single CSE query, return raw result items."""
    if not settings.google_api_key or not cx:
        return []

    await rate_limit.google_cse.acquire()
    params = {
        "key": settings.google_api_key,
        "cx": cx,
        "q": query,
        "num": 10,
    }
    try:
        resp = await client.get(_CSE_URL, params=params, timeout=20.0)
        resp.raise_for_status()
        return resp.json().get("items", [])
    except httpx.HTTPStatusError as e:
        log.warning("Google CSE HTTP error (cx=%s): %s", cx[:12], e)
        return []
    except Exception as e:
        log.warning("Google CSE error: %s", e)
        return []


class GoogleCSESource(PMIDSource):
    source = Source.GOOGLE

    async def query(self, lvg: LVGResult) -> SourceResult:
        if not settings.google_api_key:
            return SourceResult(
                source=self.source, pmids=[], error="GOOGLE_API_KEY not configured"
            )

        cache_key = f"google:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            return SourceResult(source=self.source, pmids=cached, cached=True)

        query = _build_query(lvg)
        log.debug("Google CSE query: %r", query)

        all_pmids: set[int] = set()
        async with httpx.AsyncClient() as client:
            for cx in (settings.google_cse_whitelist_id, settings.google_cse_schema_id):
                if not cx:
                    continue
                items = await _cse_query(query, cx, client)
                all_pmids.update(_extract_pmids_from_results(items))

        result = sorted(all_pmids)
        await cache_set(cache_key, result, ttl=settings.cache_ttl_google)
        return SourceResult(source=self.source, pmids=result)
