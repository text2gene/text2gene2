"""
Europe PMC source — open-access literature including preprints.

REST API: https://europepmc.org/RestfulWebService
No API key required. Covers PubMed + PMC + preprints + gray literature.
Especially useful for recent papers not yet indexed in PubMed.

Query strategy: use the expanded query builder from pipeline/expand.py which
generates gene synonym + multi-form variant OR-queries. Falls back to raw
input HGVS if expansion produces nothing.
"""
import logging

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2.pipeline.expand import build_europepmc_query
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


class EuropePMCSource(PMIDSource):
    source = Source.EUROPEPMC

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"europepmc:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            return SourceResult(source=self.source, pmids=cached, cached=True)

        query = await build_europepmc_query(lvg)
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
        return SourceResult(source=self.source, pmids=result, query_used=query)
