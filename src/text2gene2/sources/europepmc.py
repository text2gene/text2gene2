"""
Europe PMC source — open-access literature including preprints.

REST API: https://europepmc.org/RestfulWebService
No API key required. Covers PubMed + PMC + preprints + gray literature.
Especially useful for recent papers not yet indexed in PubMed.

We search with the gene symbol + protein change, returning PubMed-linked PMIDs.
"""
import logging

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def _build_query(lvg: LVGResult) -> str:
    parts = []
    if lvg.gene_symbol:
        parts.append(f'"{lvg.gene_symbol}"')

    # Use a compact protein form as a keyword — strip prefix and complex suffixes
    # so "NP_009225.1:p.(Q1756Pfs*74)" → "Q1756" (position-level, avoids wildcard issues)
    for p in lvg.hgvs_p[:1]:
        short = p.split(":")[-1] if ":" in p else p
        # Remove p.( prefix variants
        short = short.lstrip("p.(")
        # Cut off at first non-AA character sequence (frameshift marker, asterisk, etc.)
        import re as _re
        m = _re.match(r"([A-Z][a-z]{0,2}\d+)", short)
        if m:
            parts.append(f'"{m.group(1)}"')
        else:
            parts.append(f'"{short}"')

    if not parts:
        parts.append(f'"{lvg.input_hgvs}"')

    return " AND ".join(parts)


class EuropePMCSource(PMIDSource):
    source = Source.EUROPEPMC

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"europepmc:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            return SourceResult(source=self.source, pmids=cached, cached=True)

        query = _build_query(lvg)
        log.debug("Europe PMC query: %r", query)

        await rate_limit.europepmc.acquire()
        params = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": 50,
            "sort": "CITED desc",
        }

        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(_BASE, params=params, timeout=20.0)
                resp.raise_for_status()
                data = resp.json()
        except Exception as e:
            log.warning("Europe PMC error: %s", e)
            return SourceResult(source=self.source, pmids=[], error=str(e))

        pmids: list[int] = []
        for hit in data.get("resultList", {}).get("result", []):
            # Only include results with a real PubMed ID
            pmid = hit.get("pmid")
            if pmid:
                try:
                    pmids.append(int(pmid))
                except ValueError:
                    pass

        result = sorted(set(pmids))
        await cache_set(cache_key, result, ttl=settings.cache_ttl_europepmc)
        return SourceResult(source=self.source, pmids=result)
