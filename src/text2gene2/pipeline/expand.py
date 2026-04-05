"""
Query expansion for variant search — generates alternate notation forms and gene synonyms.

Ported from the original text2gene's googlequery.py + metavariant's VariantComponents.
The key insight: literature uses many notations for the same variant. We generate all
plausible forms so text-search sources (EuropePMC, Google CSE) can find papers that
use non-standard notation.

Three expansion layers:
  1. HGVS short forms — strip transcript prefixes to get bare c./p. notation
  2. Posedit slang — alternate notation for the same edit (C891T vs 891C>T vs 891C->T)
  3. Gene synonyms — aliases and previous symbols from HGNC

Usage:
    from text2gene2.pipeline.expand import expand_variant, get_gene_synonyms

    forms = expand_variant(lvg)
    # forms.coding_short    → ["c.891C>T"]
    # forms.protein_short   → ["p.Arg297Cys"]
    # forms.protein_1letter → ["R297C"]
    # forms.slang           → ["891C->T", "891C-->T", "891C/T", "C891T"]
    # forms.rsids           → ["rs12345"]

    synonyms = await get_gene_synonyms("ACVRL1")
    # ["ACVRL1", "ALK1", "HHT2", "HHT", "ACVRLK1", "ORW2"]
"""
import asyncio
import logging
import re

from text2gene2.models import LVGResult

log = logging.getLogger(__name__)

# ── Amino acid maps (from metavariant) ───────────────────────────────────────

AMINO_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

# Substitution arrow slang (from metavariant's official_to_slang_map)
_SUB_SLANG = ["->", "-->", "/"]

# Ter synonyms used in literature
_TER_SLANG = ["*", "X", "Stop"]

# Gene synonyms we never use — too ambiguous, pollute search results.
# Common English words, abbreviations, amino acids, or disease names (not gene-specific).
_GENE_SYNONYM_BLACKLIST = {
    # Common English words that are gene aliases
    "END", "ALL", "NOT", "SET", "REST", "CAT", "FAT", "MAT", "HIT", "KIT",
    "BAD", "BAG", "TAG", "MAP", "GAP", "CAP", "RAP", "SAP", "TAP",
    "MAD", "TIP", "SIP", "HIP", "DIP", "RIP", "ZIP", "TOP",
    "FIT", "BIT", "PIT", "HAT", "RAT", "SAT",
    # Too-short or ambiguous scientific abbreviations
    "CF",      # "cf." = "compare" in academic writing
    "MET",     # methionine amino acid
    "DEC",     # original blacklist entry
    "MIS",     # common word
    "CAMP",    # cyclic AMP
    "CHIP",    # chromatin immunoprecipitation
}

# Minimum length for gene synonyms — anything shorter than 3 chars is too ambiguous
_MIN_SYNONYM_LENGTH = 3

# Regex to parse protein changes: p.(Arg297Cys) or p.Arg297Cys
_RE_PROTEIN = re.compile(
    r"p\.?\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=|\?|fs)\)?"
)

# Regex to parse coding changes: c.891C>T, c.4964_4982del, c.6869dupA, etc.
_RE_CODING_SUB = re.compile(r"c\.(\d+(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?)([ACGT])>([ACGT])")
_RE_CODING_DEL = re.compile(r"c\.(\d+(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?)del([ACGT]*)")
_RE_CODING_DUP = re.compile(r"c\.(\d+(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?)dup([ACGT]*)")
_RE_CODING_INS = re.compile(r"c\.(\d+(?:[+-]\d+)?_\d+(?:[+-]\d+)?)ins([ACGT]+)")
_RE_CODING_DELINS = re.compile(r"c\.(\d+(?:[+-]\d+)?(?:_\d+(?:[+-]\d+)?)?)del[ACGT]*ins([ACGT]+)")


class ExpandedVariant:
    """All expanded forms of a variant, ready for query building."""

    def __init__(self):
        self.gene_symbol: str | None = None
        self.coding_short: list[str] = []      # e.g. ["c.891C>T"]
        self.protein_short: list[str] = []     # e.g. ["p.Arg297Cys", "p.(Arg297Cys)"]
        self.protein_1letter: list[str] = []   # e.g. ["R297C"]
        self.genomic_short: list[str] = []     # e.g. ["g.12345C>T"]
        self.slang: list[str] = []             # e.g. ["891C->T", "C891T"]
        self.rsids: list[str] = []             # e.g. ["rs12345"]
        self.all_hgvs: list[str] = []          # full HGVS with transcript

    @property
    def all_search_forms(self) -> list[str]:
        """All unique search terms, deduplicated, ordered by specificity."""
        seen = set()
        result = []
        for form in (
            self.coding_short + self.protein_short + self.protein_1letter +
            self.slang + self.rsids + self.genomic_short
        ):
            if form and form not in seen:
                seen.add(form)
                result.append(form)
        return result


def _strip_transcript(hgvs: str) -> str:
    """NM_000492.4:c.891C>T → c.891C>T"""
    return hgvs.split(":")[-1] if ":" in hgvs else hgvs


_RE_PROTEIN_DEL = re.compile(
    r"p\.?\(?([A-Z][a-z]{2})(\d+)del\)?"
)
_RE_PROTEIN_FS = re.compile(
    r"p\.?\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?fs"
)


def _protein_to_1letter(p_short: str) -> str | None:
    """p.Arg297Cys → R297C, p.(Phe508del) → F508del, p.(Arg297CysfsTer5) → R297fs"""
    # Frameshift first (before substitution, since p.Arg297Cysfs matches both)
    m = _RE_PROTEIN_FS.match(p_short)
    if m:
        ref1 = AMINO_3TO1.get(m.group(1))
        if ref1:
            return f"{ref1}{m.group(2)}fs"

    # Deletion: p.(Phe508del) → F508del
    m = _RE_PROTEIN_DEL.match(p_short)
    if m:
        ref1 = AMINO_3TO1.get(m.group(1))
        if ref1:
            return f"{ref1}{m.group(2)}del"

    # Standard substitution
    m = _RE_PROTEIN.match(p_short)
    if m:
        ref3, pos, alt3 = m.group(1), m.group(2), m.group(3)
        ref1 = AMINO_3TO1.get(ref3)
        alt1 = AMINO_3TO1.get(alt3, alt3)  # keep "=", "?" as-is
        if ref1 is None:
            return None
        return f"{ref1}{pos}{alt1}"

    return None


def _generate_sub_slang(pos: str, ref: str, alt: str) -> list[str]:
    """Generate alternate substitution notations.

    We keep the c. prefix on arrow forms to avoid matching genotype calls.
    The ref-pos-alt form (e.g. C891T) is common in older literature and
    specific enough with a gene AND-clause to be safe.
    """
    out = []
    # c.891C->T, c.891C-->T (keep c. prefix — bare "891C/T" looks like a genotype)
    for arrow in ("->", "-->"):
        out.append(f"c.{pos}{ref}{arrow}{alt}")
    # C891T (ref-pos-alt, common in older literature, no c. prefix is how it appears)
    out.append(f"{ref}{pos}{alt}")
    return out


def _generate_protein_slang(p_short: str) -> list[str]:
    """Generate alternate protein notations."""
    out = []
    # Strip parentheses for matching: p.(Arg297Cys) → p.Arg297Cys
    stripped = p_short.replace("(", "").replace(")", "")
    if stripped != p_short:
        out.append(stripped)

    # Ter → *, X, Stop
    for ter_form in _TER_SLANG:
        replaced = stripped.replace("Ter", ter_form)
        if replaced != stripped:
            out.append(replaced)

    # Truncate frameshift: p.Arg297GlyfsTer5 → p.Arg297fs
    fs_pos = stripped.find("fs")
    if fs_pos > -1:
        out.append(stripped[:fs_pos + 2])

    return out


def _generate_coding_slang(c_short: str) -> list[str]:
    """Generate alternate coding notations."""
    out = []

    # Substitution: c.891C>T → 891C->T, 891C-->T, 891C/T, C891T
    m = _RE_CODING_SUB.match(c_short)
    if m:
        out.extend(_generate_sub_slang(m.group(1), m.group(2), m.group(3)))
        return out

    # Deletion with bases: c.4964_4982delCTGGCC → 4964_4982del (short form)
    m = _RE_CODING_DEL.match(c_short)
    if m and m.group(2):  # has bases specified
        out.append(f"c.{m.group(1)}del")
        return out

    # Duplication with bases: c.6869dupA → 6869dup
    m = _RE_CODING_DUP.match(c_short)
    if m and m.group(2):
        out.append(f"c.{m.group(1)}dup")
        return out

    return out


def expand_variant(lvg: LVGResult) -> ExpandedVariant:
    """
    Generate all expanded search forms from an LVGResult.

    Takes the VariantValidator output and produces:
    - Short forms (transcript stripped)
    - 1-letter protein codes
    - Posedit slang (arrow variants, ref-pos-alt, short deletions)
    - Protein slang (Ter variants, frameshift truncation)
    """
    exp = ExpandedVariant()
    exp.gene_symbol = lvg.gene_symbol
    exp.rsids = list(lvg.rsids)
    exp.all_hgvs = lvg.all_hgvs

    seen_coding = set()
    seen_protein = set()

    # Coding forms
    for h in lvg.hgvs_c:
        short = _strip_transcript(h)
        if short not in seen_coding:
            seen_coding.add(short)
            exp.coding_short.append(short)
            exp.slang.extend(_generate_coding_slang(short))

    # Protein forms
    for h in lvg.hgvs_p:
        short = _strip_transcript(h)
        if short not in seen_protein:
            seen_protein.add(short)
            exp.protein_short.append(short)
            exp.slang.extend(_generate_protein_slang(short))

            # 1-letter code
            one_letter = _protein_to_1letter(short)
            if one_letter and one_letter not in exp.protein_1letter:
                exp.protein_1letter.append(one_letter)

    # Genomic forms (less useful for text search, but include)
    seen_genomic = set()
    for h in lvg.hgvs_g:
        short = _strip_transcript(h)
        if short not in seen_genomic:
            seen_genomic.add(short)
            exp.genomic_short.append(short)

    # Deduplicate slang
    seen_slang = set()
    deduped = []
    for s in exp.slang:
        if s not in seen_slang:
            seen_slang.add(s)
            deduped.append(s)
    exp.slang = deduped

    return exp


# ── Gene synonyms from HGNC ───���─────────────────────────────────────────────

_synonym_cache: dict[str, list[str]] = {}


def _fetch_gene_synonyms_sync(symbol: str) -> list[str]:
    """Query hgnc.gene for alias_symbols and prev_symbols."""
    from text2gene2.db import get_medgen_conn, reset_medgen_conn

    conn = get_medgen_conn()
    if conn is None:
        return [symbol]

    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT alias_symbols, prev_symbols FROM hgnc.gene WHERE symbol = %s",
                (symbol,)
            )
            row = cur.fetchone()
    except Exception as e:
        log.warning("expand: HGNC query error for %s: %s", symbol, e)
        reset_medgen_conn()
        return [symbol]

    if not row:
        return [symbol]

    synonyms = [symbol]
    for field in row:
        if field:
            for alias in field.split("|"):
                alias = alias.strip()
                if (alias
                    and alias not in _GENE_SYNONYM_BLACKLIST
                    and len(alias) >= _MIN_SYNONYM_LENGTH
                    and alias not in synonyms):
                    synonyms.append(alias)
    return synonyms


async def get_gene_synonyms(symbol: str) -> list[str]:
    """
    Return gene symbol + all known aliases/previous symbols from HGNC.
    Cached in memory (gene names rarely change).
    """
    if symbol in _synonym_cache:
        return _synonym_cache[symbol]

    synonyms = await asyncio.to_thread(_fetch_gene_synonyms_sync, symbol)
    _synonym_cache[symbol] = synonyms
    return synonyms


# ── Query builders for text-search sources ───────────────────────────────────

async def build_europepmc_query(lvg: LVGResult) -> str:
    """
    Build an expanded EuropePMC query from an LVGResult.

    Strategy: (gene OR alias1 OR alias2) AND (variant_form1 OR variant_form2 OR ...)

    EuropePMC supports boolean queries with OR/AND and quoted phrases.
    """
    exp = expand_variant(lvg)

    if not exp.gene_symbol:
        # No gene — fall back to quoted input HGVS
        return f'"{lvg.input_hgvs}"'

    # Gene clause with synonyms
    synonyms = await get_gene_synonyms(exp.gene_symbol)
    if len(synonyms) > 1:
        gene_clause = "(" + " OR ".join(f'"{s}"' for s in synonyms) + ")"
    else:
        gene_clause = f'"{synonyms[0]}"'

    # Variant clause — pick the most specific forms, avoid bloating the query
    variant_terms = []

    # Protein forms are most commonly found in text
    for p in exp.protein_short[:2]:
        variant_terms.append(f'"{p}"')
    for p1 in exp.protein_1letter[:2]:
        variant_terms.append(f'"{p1}"')

    # Coding short forms
    for c in exp.coding_short[:2]:
        variant_terms.append(f'"{c}"')

    # Slang forms (most valuable — catches older notation)
    for s in exp.slang[:4]:
        variant_terms.append(f'"{s}"')

    # rsIDs
    for rs in exp.rsids[:2]:
        variant_terms.append(f'"{rs}"')

    if not variant_terms:
        return f'"{lvg.input_hgvs}"'

    variant_clause = "(" + " OR ".join(variant_terms) + ")"
    return f"{gene_clause} AND {variant_clause}"
