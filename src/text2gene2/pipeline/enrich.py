"""
Article metadata enrichment for CitationTable.

Batch-fetches PubMed metadata (title, authors, journal, year, DOI, PMC ID)
for all citations via NCBI efetch and populates the Citation model fields.

Results are cached in Redis per PMID with a long TTL (metadata rarely changes).
This is run automatically for web UI requests but skipped for the JSON API
(use ?enrich=true to opt in).
"""
import json
import logging
from xml.etree import ElementTree

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.models import Citation, CitationTable
from text2gene2 import rate_limit
from text2gene2.config import settings

log = logging.getLogger(__name__)

_EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_BATCH   = 100
_TTL     = 60 * 60 * 24 * 30   # 30 days — article metadata is stable


def _api_params() -> dict:
    p: dict = {}
    if settings.ncbi_api_key:
        p["api_key"] = settings.ncbi_api_key
    return p


def _first_author(article) -> str:
    authors = article.findall(".//Author")
    if not authors:
        return ""
    a = authors[0]
    last  = (a.findtext("LastName") or "").strip()
    initials = (a.findtext("Initials") or "").strip()
    first = (a.findtext("ForeName") or initials).strip()
    name = f"{last} {first}".strip() if first else last
    return f"{name} et al." if len(authors) > 1 else name


def _parse_article(article) -> dict:
    """Extract metadata fields from a PubmedArticle XML element."""
    meta: dict = {}

    # Title — strip trailing period, handle MathML/italics fragments
    title_el = article.find(".//ArticleTitle")
    if title_el is not None:
        meta["title"] = "".join(title_el.itertext()).strip().rstrip(".")

    # Authors
    meta["authors"] = _first_author(article)

    # Journal abbreviation
    j = article.find(".//Journal")
    if j is not None:
        meta["journal"] = (
            j.findtext("ISOAbbreviation")
            or j.findtext("Title")
            or ""
        ).strip()

    # Year — prefer MedlineDate → Year → PubDate/MedlineDate
    for xpath in [".//PubDate/Year", ".//PubMedPubDate[@PubStatus='pubmed']/Year"]:
        yr = article.findtext(xpath)
        if yr and yr.isdigit():
            meta["year"] = int(yr)
            break

    # DOI and PMC
    for id_el in article.findall(".//ArticleId"):
        id_type = id_el.get("IdType", "")
        val = (id_el.text or "").strip()
        if id_type == "doi" and val:
            meta["doi"] = val
        elif id_type == "pmc" and val:
            meta["pmc"] = val

    # Abstract (first 400 chars — for display)
    parts = [el.text for el in article.findall(".//AbstractText") if el.text]
    if parts:
        full = " ".join(parts)
        meta["abstract_snippet"] = full[:400] + ("…" if len(full) > 400 else "")

    return meta


async def _fetch_meta_batch(pmids: list[int]) -> dict[int, dict]:
    """Fetch metadata for a batch of PMIDs. Returns {pmid: meta_dict}."""
    result: dict[int, dict] = {}
    params = {
        **_api_params(),
        "db": "pubmed",
        "id": ",".join(str(p) for p in pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    await rate_limit.clinvar.acquire()
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(f"{_EUTILS}/efetch.fcgi", params=params, timeout=30.0)
            resp.raise_for_status()
            tree = ElementTree.fromstring(resp.text)
            for article in tree.findall(".//PubmedArticle"):
                pmid_el = article.find(".//PMID")
                if pmid_el is None or not pmid_el.text:
                    continue
                pmid = int(pmid_el.text)
                result[pmid] = _parse_article(article)
    except Exception as e:
        log.warning("enrich efetch error for batch %d…: %s", pmids[0], e)
    return result


async def enrich_citations(table: CitationTable) -> CitationTable:
    """
    Populate Citation metadata fields (title, authors, journal, year, doi)
    for all citations in the table. Uses Redis cache; only fetches uncached PMIDs.
    """
    if not table.citations:
        return table

    # Split into cached / uncached
    meta_map: dict[int, dict] = {}
    uncached: list[int] = []

    for c in table.citations:
        key = f"meta:{c.pmid}"
        cached = await cache_get(key)
        if cached is not None:
            meta_map[c.pmid] = cached
        else:
            uncached.append(c.pmid)

    # Batch-fetch uncached
    for i in range(0, len(uncached), _BATCH):
        batch = uncached[i : i + _BATCH]
        fetched = await _fetch_meta_batch(batch)
        for pmid, meta in fetched.items():
            meta_map[pmid] = meta
            await cache_set(f"meta:{pmid}", meta, ttl=_TTL)
        # Mark misses so we don't re-fetch
        for pmid in batch:
            if pmid not in fetched:
                await cache_set(f"meta:{pmid}", {}, ttl=_TTL)

    # Apply to citations
    for citation in table.citations:
        meta = meta_map.get(citation.pmid, {})
        if meta.get("title"):
            citation.title = meta["title"]
        if meta.get("authors"):
            citation.authors = meta["authors"]
        if meta.get("journal"):
            citation.journal = meta["journal"]
        if meta.get("year"):
            citation.year = meta["year"]
        if meta.get("doi"):
            citation.doi = meta["doi"]
        # pmc and abstract_snippet go into extras for the template
        citation.pmc      = meta.get("pmc")
        citation.abstract_snippet = meta.get("abstract_snippet")

    return table
