"""
Europe PMC source — open-access literature including preprints.

REST API: https://europepmc.org/RestfulWebService
No API key required. Covers PubMed + PMC + preprints + gray literature.
Especially useful for recent papers not yet indexed in PubMed.

Query strategy (in priority order):
  1. gene + protein change (when protein effect is meaningful, e.g. p.Arg408Trp)
  2. gene + c. short form (for splice-site, frameshifts, and p.? variants)
  3. raw input HGVS (last resort)

Returning just the gene name is never acceptable — it produces thousands of
off-topic results (e.g. searching "CFTR" returns Hsp70 chaperone papers).
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

_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Protein changes that carry no specific information about the variant
_UNINFORMATIVE_P = {"?", "=", "0", "0?"}


def _protein_keyword(hgvs_p_list: list[str]) -> str | None:
    """
    Extract a meaningful protein-level search keyword from the first hgvs_p entry.
    Returns None if the protein effect is unknown/uninformative (p.?, p.=, splice).
    """
    for p in hgvs_p_list[:1]:
        short = p.split(":")[-1] if ":" in p else p
        short = short.lstrip("p.(").rstrip(")")
        if short in _UNINFORMATIVE_P or not short:
            return None
        # Extract amino-acid + position (e.g. "Arg408" from "Arg408Trp")
        m = re.match(r"([A-Z][a-z]{0,2}\d+)", short)
        return f'"{m.group(1)}"' if m else None
    return None


def _coding_keyword(hgvs_c_list: list[str]) -> str | None:
    """
    Extract the short c. form (without transcript) as a search keyword.
    e.g. "NM_000492.4:c.3963+1G>A" → '"c.3963+1G>A"'
    """
    for h in hgvs_c_list[:1]:
        short = h.split(":")[-1] if ":" in h else h
        if short.startswith("c.") and len(short) > 2:
            return f'"{short}"'
    return None


def _build_query(lvg: LVGResult) -> str:
    gene = f'"{lvg.gene_symbol}"' if lvg.gene_symbol else None

    # Priority 1: gene + meaningful protein change
    p_kw = _protein_keyword(lvg.hgvs_p)
    if gene and p_kw:
        return f"{gene} AND {p_kw}"

    # Priority 2: gene + c. notation (splice sites, frameshifts, intronic, p.?)
    c_kw = _coding_keyword(lvg.hgvs_c)
    if gene and c_kw:
        return f"{gene} AND {c_kw}"

    # Priority 3: bare input HGVS — at least variant-specific
    return f'"{lvg.input_hgvs}"'


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
