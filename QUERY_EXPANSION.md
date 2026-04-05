# Query Expansion in text2gene2

## Overview

Query expansion generates multiple search forms for a single genetic variant so
text-search sources (EuropePMC, Google CSE) can find papers that use non-standard
notation. Literature uses many ways to describe the same variant — HGVS official
notation, shorthand, 1-letter amino acid codes, old arrow syntax, gene aliases —
and we need to search for all of them.

The implementation lives in `src/text2gene2/pipeline/expand.py`.

## Provenance

Ported from the original text2gene project, which had a mature query expansion
system across three components:

- **`metavariant/metavariant/components.py`** — `VariantComponents` class with
  `posedit_slang` property generating alternate notation for substitutions,
  deletions, duplications, frameshifts, and protein changes. This was the core
  expansion engine.

- **`text2gene/text2gene/googlequery.py`** — `GoogleCSEngine.build_query()` which
  combined gene synonyms (via `medgen.annotate.gene.GeneSynonyms`) with posedit
  slang into OR-joined queries, with term counting against Google's 32-term limit.

- **`text2gene/text2gene/pmid_lookups.py`** — Used LVG expansion across all
  sequence types (c, p, g, n) for both ClinVar and PubTator lookups.

text2gene2 reimplements the expansion logic directly (no metavariant dependency)
and pulls gene synonyms from the local `hgnc.gene` table on loki Postgres instead
of the old `medgen.annotate.gene.GeneSynonyms` library.

## Three Expansion Layers

### 1. HGVS Short Forms

Strip transcript/protein accession prefixes to get bare notation:

| Input (from VariantValidator) | Short form |
|-------------------------------|------------|
| `NM_000020.3:c.1120C>T` | `c.1120C>T` |
| `NP_000011.2:p.(Arg374Trp)` | `p.(Arg374Trp)` |
| `NC_000012.12:g.51916568G>A` | `g.51916568G>A` |

These are exact HGVS notation without the accession. Safe and unambiguous.

### 2. Posedit Slang

Generate alternate notations for the same edit, based on patterns found in
published literature. Each edit type has its own slang rules:

#### Substitutions (`c.1120C>T`)

| Slang form | Notes |
|-----------|-------|
| `c.1120C->T` | Arrow notation (common in older papers) |
| `c.1120C-->T` | Double-arrow variant |
| `C1120T` | Ref-pos-alt (no prefix, compact form used pre-HGVS standardization) |

**Dropped:** `1120C/T` — looks like a SNP genotype call (e.g., from HapMap/dbSNP
genotyping data). The `/` notation is ambiguous: `C/T` at position 1120 could mean
heterozygous genotype, not a specific mutation. Removed to avoid false positives
from genotyping papers.

**Decision:** Arrow slang retains the `c.` prefix (`c.1120C->T`, not `1120C->T`)
because bare numeric forms match too broadly in full-text search (figure labels,
table entries, clone IDs). The ref-pos-alt form (`C1120T`) omits `c.` because
that's genuinely how it appears in older papers — and the gene AND-clause provides
sufficient constraint.

#### Deletions (`c.1521_1523delCTT`)

| Slang form | Notes |
|-----------|-------|
| `c.1521_1523del` | Short form — bases omitted (common shorthand) |

The long form with explicit bases is rarely written out in text. Most papers just
say "del" without specifying the deleted sequence.

#### Duplications (`c.6869dupA`)

| Slang form | Notes |
|-----------|-------|
| `c.6869dup` | Short form — duplicated base(s) omitted |

Same logic as deletions.

#### Insertions and Indels

No slang forms generated. These are complex enough that the official notation is
what appears in papers. Could revisit if we find systematic alternate patterns.

#### Protein Changes (`p.(Arg374Trp)`)

| Slang form | Notes |
|-----------|-------|
| `p.Arg374Trp` | Parentheses stripped (some papers omit them) |
| `R374W` | 1-letter amino acid code (very common in literature) |

For stop codons, Ter is expanded to synonyms:

| Input | Slang forms |
|-------|------------|
| `p.(Gly12Ter)` | `p.Gly12*`, `p.Gly12X`, `p.Gly12Stop` |

For frameshifts, downstream notation is truncated:

| Input | Slang forms |
|-------|------------|
| `p.(Glu1309AspfsTer4)` | `p.Glu1309Aspfs` (truncated), plus Ter→`*`/`X`/`Stop` |

1-letter protein codes are also generated for deletions and frameshifts:

| Input | 1-letter |
|-------|----------|
| `p.(Phe508del)` | `F508del` |
| `p.(Arg297CysfsTer5)` | `R297fs` |

#### Genomic Coordinates

Genomic HGVS forms (`g.51916568G>A`) are included as short forms but no slang
is generated. The original text2gene experimentally determined that genomic
coordinates do not help in text search — they are rarely cited in papers. We
include them for completeness but don't prioritize them.

### 3. Gene Synonyms

Pulled from `hgnc.gene` on loki Postgres (`alias_symbols` and `prev_symbols`
columns, pipe-delimited). Used to build OR-joined gene clauses so papers using
old or alternate gene names are found.

Example: `ACVRL1` → `ACVRL1`, `HHT2`, `ALK1`, `HHT`, `ACVRLK1`, `ORW2`

#### Blacklist

Some gene aliases are common English words, scientific abbreviations, or disease
names that would match far too broadly in full-text search. These are blacklisted:

| Blacklisted | Gene | Reason |
|-------------|------|--------|
| `END` | ENG | Common English word |
| `CF` | CFTR | "cf." means "compare" in academic writing |
| `MET` | MET | Methionine (amino acid), appears in every proteomics paper |
| `ALL` | BCR | Common English word |
| `NOT`, `SET`, `REST` | various | Common English words |
| `CAT`, `FAT`, `MAT`, `HIT` | various | Common English words |
| `BAD`, `BAG`, `TAG`, `MAP` | various | Common English / abbreviations |
| `GAP`, `CAP`, `RAP`, `SAP`, `TAP` | various | Common abbreviations |
| `MAD`, `TIP`, `SIP`, `HIP`, `DIP`, `RIP`, `ZIP` | various | Common words |
| `TOP`, `FIT`, `BIT`, `PIT`, `HAT`, `RAT`, `SAT` | various | Common words |
| `MIS` | AMH | Common word |
| `DEC` | DEC1 | Original text2gene blacklist |
| `CAMP` | CHAMP1 | Cyclic AMP — ubiquitous in biology papers |
| `CHIP` | STUB1 | Chromatin immunoprecipitation — ubiquitous method |

#### Minimum Length

All synonyms shorter than 3 characters are rejected regardless of blacklist.
Two-letter gene aliases (e.g., `DP` for some APC aliases) are too ambiguous.

#### Synonyms We Keep Despite Being Short

The following are 3-character gene aliases that pass the filter. They're the
**actual approved HGNC symbols**, so we can't drop them:

- `ENG`, `APC`, `MET`, `KIT`, `RET` — primary gene symbols, not aliases
- `HHT` — disease abbreviation but specific to hereditary hemorrhagic telangiectasia;
  acceptable risk when AND-joined with variant notation
- `ORW` — Osler-Rendu-Weber, disease name; similar reasoning
- `DP2`, `DP3`, `PBT`, `PTC` — obscure enough to not collide

The AND-conjunction with specific variant notation is the main safeguard: you'd
need a paper that mentions both "HHT" AND "c.1120C->T" to get a false positive.

## How Queries Are Built

### EuropePMC (`build_europepmc_query`)

Produces a boolean query:

```
("gene1" OR "gene2" OR ...) AND ("variant_form1" OR "variant_form2" OR ...)
```

Variant forms are selected in priority order (most specific first), capped to
avoid query bloat:

1. Protein short forms (up to 2)
2. 1-letter protein codes (up to 2)
3. Coding short forms (up to 2)
4. Slang forms (up to 4)
5. rsIDs (up to 2)

Total variant terms: up to ~14 per query. Gene synonyms: unlimited (typically 3-8).

### Google CSE (not yet wired — see issue #17)

Will need a different builder due to Google's 32-term limit. The original text2gene
counted terms carefully (Google splits on special characters, so `891C>T` counts
as 2 terms). When implementing, reference `text2gene/text2gene/googlequery.py`
for the term-counting logic.

### LOVD Matching (planned)

The expanded forms will also be used to match variants against LOVD database
entries, which store variants in their own notation. The `ExpandedVariant` object
is source-agnostic — any consumer can use `all_search_forms` or pick specific
form lists.

## Design Philosophy: Precision Over Recall

The goal of text2gene2 is to produce a list of literature that a clinician or
curator can **trust without manually re-verifying every result**. A false positive
— a paper that doesn't actually discuss the queried variant — is worse than a
missed paper, because it wastes expert review time and can mislead clinical
interpretation. We are not building a search engine where users expect to sift
through noisy results. We are building an evidence aggregator whose output should
be actionable.

This doesn't mean we avoid broad queries — it means we **rank and annotate**
results so users can tell the difference between a high-confidence hit and a
speculative one. It's fine to return "maybe" results from looser expansion
strategies, as long as:

1. Every result carries metadata about **how it was found** (which source, which
   query form matched, exact vs. slang vs. synonym)
2. Results are **sorted by confidence** — a PMID confirmed by 3 independent
   sources and validated against the abstract ranks above a PMID found only
   by an LLM-expanded disease query
3. The UI clearly communicates **why a result might be shaky** — e.g., "found via
   gene synonym 'ALK1' + 1-letter protein code 'R374W', not confirmed in abstract"

The goal is not to hide uncertain results but to present them honestly so the
user can make informed decisions about what to review. A clinician glancing at
search results should immediately see which papers are solid and which need a
closer look. The expansion system's job is to cast a wide net; the validation
and ranking systems' job is to sort the catch.

This matters most for exactly the cases text2gene2 is built for: **rare diseases
and difficult curation cases**. A well-studied variant in BRCA1 will have dozens
of easily-found papers — you don't need query expansion for that. But a novel
missense variant in ACVRL1 with zero LitVar2 hits might have one relevant paper
from 2004 that used the old gene name "ALK-1" and described the mutation as
"C1120T" without HGVS notation. That paper could be the difference between a VUS
and a likely pathogenic classification. The entire point of expansion is to find
evidence that standard searches miss, and the entire point of confidence ranking
is to surface those finds without burying the user in noise.

## False Positive Analysis

Every expansion form is a tradeoff between recall (finding more papers) and
precision (avoiding irrelevant matches). Here's the risk assessment:

| Form | Example | False Positive Risk | Why |
|------|---------|-------------------|-----|
| Full HGVS c./p. | `c.1120C>T` | None | Unambiguous standard notation |
| rsID | `rs121909660` | None | Globally unique identifier |
| Arrow slang (with c.) | `c.1120C->T` | Very low | Prefix constrains context |
| Ref-pos-alt | `C1120T` | Low | Rare pattern outside genetics |
| 1-letter protein | `R374W` | Low | Rare pattern in non-genetics text |
| Deletion short | `c.1521_1523del` | Very low | Has c. prefix |
| Gene synonyms | `ALK1`, `HHT2` | Low | Blacklisted common words |

**Key safeguard:** All forms are used in AND-conjunction with gene names. A false
positive requires matching both gene AND variant notation in the same paper.

## What We Explicitly Don't Do

- **No fuzzy matching.** All expansions are deterministic string transforms. We
  never do edit-distance matching, phonetic matching, or probabilistic scoring.

- **No LLM expansion (yet).** Issue #15 tracks LLM-assisted query expansion for
  zero-result variants. When implemented, LLM-generated terms should be treated
  as lower-confidence and validated more strictly.

- **No genomic coordinate expansion.** Experimentally proven unhelpful in the
  original text2gene. Papers almost never cite `g.` coordinates.

- **No bare numeric slang.** Forms like `1120C>T` (without `c.` prefix) are not
  generated because they match too broadly.

- **No gene-only queries.** If expansion produces no variant terms, we fall back
  to the raw input HGVS, never to just the gene name. A gene-name-only query
  returns thousands of irrelevant results.
