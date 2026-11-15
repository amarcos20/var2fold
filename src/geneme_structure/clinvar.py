"""Find known clinical variants for a gene in ClinVar (NCBI E-utilities).

`ClinVarClient.variants_for_gene("BRCA1")` returns GRCh38 variants ready for the
pipeline, so a user can start from a gene symbol instead of coordinates.

* Search: esearch on db=clinvar, filtered by significance, consequence and
  `single_gene` (variants spanning several genes are excluded).
* Coordinates: each record's canonical SPDI (GRCh38 RefSeq, 0-based) is
  converted to a VCF-style 1-based variant. Pure insertions/deletions need a
  left anchor base, fetched from the Ensembl sequence endpoint.
* Ranking: review status (stars) first, then Pathogenic before Likely
  pathogenic. Ties at the cut-off are sampled evenly along the gene, so a
  small selection still covers the whole protein instead of one end of it.
"""

from __future__ import annotations

import logging
import re
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import requests

from geneme_structure import __version__
from geneme_structure.http_session import build_session
from geneme_structure.models import GenomicVariant, normalize_chromosome

logger = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_REST = "https://rest.ensembl.org"
ESUMMARY_BATCH = 200
MAX_SEARCH_IDS = 5_000

# ClinVar "molecular consequence" terms, grouped into user-facing sets.
CONSEQUENCE_SETS: dict[str, tuple[str, ...]] = {
    "missense": ("missense variant",),
    "coding": (
        "missense variant",
        "nonsense",
        "frameshift variant",
        "inframe deletion",
        "inframe insertion",
        "inframe indel",
    ),
    "splice": ("splice donor variant", "splice acceptor variant"),
    "any": (),
}
SIGNIFICANCE_FILTERS = {
    "pathogenic": "clinsig_pathogenic",
    "likely_pathogenic": "clinsig_likely_pathogenic",
}
# One amino-acid substitution, e.g. R175H (excludes stop gains 'E1817*' and frameshifts).
MISSENSE_CHANGE = re.compile(r"^[A-Z]\d+[A-Z]$")
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
}


@dataclass(frozen=True)
class ClinVarVariant:
    variant: GenomicVariant
    accession: str  # VCV accession
    title: str
    protein_change: str | None
    classification: str
    review_status: str
    consequences: tuple[str, ...] = field(default_factory=tuple)

    @property
    def stars(self) -> int:
        return REVIEW_STARS.get(self.review_status, 0)

    @property
    def label(self) -> str:
        """Short display name: protein change (e.g. R175H) or the VCV accession."""
        return self.protein_change or self.accession


# --------------------------------------------------------------------------- pure helpers


def build_query(gene: str, consequence: str, significance: Sequence[str]) -> str:
    if consequence not in CONSEQUENCE_SETS:
        raise ValueError(f"consequence must be one of {sorted(CONSEQUENCE_SETS)}")
    unknown = set(significance) - set(SIGNIFICANCE_FILTERS)
    if unknown or not significance:
        raise ValueError(f"significance must be a subset of {sorted(SIGNIFICANCE_FILTERS)}")
    if not gene.replace("-", "").replace(".", "").isalnum():
        raise ValueError(f"{gene!r} does not look like a gene symbol")

    sig = " OR ".join(f"{SIGNIFICANCE_FILTERS[s]}[prop]" for s in significance)
    parts = [f"{gene}[gene]", f"({sig})", "single_gene[prop]"]
    terms = CONSEQUENCE_SETS[consequence]
    if terms:
        parts.append("(" + " OR ".join(f'"{t}"[molecular consequence]' for t in terms) + ")")
    return " AND ".join(parts)


def spdi_to_vcf(spdi: str) -> tuple[int, str, str]:
    """Canonical SPDI -> (1-based position, REF, ALT), trimmed VCF-style.

    REF or ALT is empty for a pure insertion/deletion: the caller must then
    prepend the reference base at position-1 (the VCF anchor).
    """
    _seq, pos0, ref, alt = spdi.split(":")
    pos = int(pos0) + 1
    if ref and alt:
        # Canonical SPDI expands repeats ("AAAA:AAA"); trim shared bases but keep >=1 each.
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _grch38_chromosome(doc: dict[str, Any]) -> str | None:
    for vs in doc.get("variation_set") or []:
        for loc in vs.get("variation_loc") or []:
            if loc.get("assembly_name") == "GRCh38" and loc.get("chr"):
                return normalize_chromosome(str(loc["chr"]))
    return None


def parse_summary(doc: dict[str, Any], gene: str) -> dict[str, Any] | None:
    """Fields needed to build a `ClinVarVariant`, or None if the record is unusable."""
    sets = doc.get("variation_set") or []
    if len(sets) != 1:  # haplotypes / compound records
        return None
    spdi = sets[0].get("canonical_spdi")
    chrom = _grch38_chromosome(doc)
    if not spdi or not chrom:
        return None
    if gene.upper() not in {g.get("symbol", "").upper() for g in doc.get("genes") or []}:
        return None
    pos, ref, alt = spdi_to_vcf(spdi)
    germline = doc.get("germline_classification") or {}
    protein = (doc.get("protein_change") or "").split(",")[0].strip() or None
    return {
        "chromosome": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "accession": doc.get("accession", ""),
        "title": doc.get("title", ""),
        "protein_change": protein,
        "classification": germline.get("description", ""),
        "review_status": germline.get("review_status", ""),
        "consequences": tuple(doc.get("molecular_consequence_list") or ()),
    }


def _tier(p: dict[str, Any]) -> tuple[int, int]:
    """Lower is better: review stars, then Pathogenic before Likely pathogenic."""
    return (
        -REVIEW_STARS.get(p["review_status"], 0),
        0 if p["classification"].lower().startswith("pathogenic") else 1,
    )


def matches_consequence(p: dict[str, Any], consequence: str) -> bool:
    """ClinVar lists consequences over all transcripts; for 'missense' also require
    the reported protein change to be a single substitution."""
    if consequence != "missense":
        return True
    change = p.get("protein_change")
    return change is None or bool(MISSENSE_CHANGE.match(change))


def select_spread(parsed: Sequence[dict[str, Any]], n: int | None) -> list[dict[str, Any]]:
    """Best tiers first; inside the tier that crosses the cut-off, pick evenly by position.

    Deterministic, and the result is returned in genomic order within each tier.
    """
    tiers: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for p in parsed:
        tiers.setdefault(_tier(p), []).append(p)
    out: list[dict[str, Any]] = []
    for key in sorted(tiers):
        members = sorted(tiers[key], key=lambda p: (p["chromosome"], p["position"]))
        room = len(members) if n is None else n - len(out)
        if room <= 0:
            break
        if room >= len(members):
            out.extend(members)
        else:
            step = len(members) / room
            out.extend(members[int(i * step + step / 2)] for i in range(room))
    return out


# --------------------------------------------------------------------------- client


class ClinVarClient:
    def __init__(
        self,
        *,
        timeout: float = 30.0,
        session: requests.Session | None = None,
        api_key: str | None = None,
    ) -> None:
        self.timeout = timeout
        self.session = session or build_session()
        self.api_key = api_key
        # NCBI asks for <=3 requests/s without a key, <=10 with one.
        self.min_interval_s = 0.11 if api_key else 0.34
        self._last_call = 0.0

    def __repr__(self) -> str:
        return f"ClinVarClient(api_key={'***' if self.api_key else None})"

    def _eutils(self, endpoint: str, params: dict[str, Any]) -> dict[str, Any]:
        wait = self.min_interval_s - (time.monotonic() - self._last_call)
        if wait > 0:
            time.sleep(wait)
        params = {"db": "clinvar", "retmode": "json", "tool": "geneme-structure", **params}
        if self.api_key:
            params["api_key"] = self.api_key
        # POST keeps long id lists out of the URL (and the key out of server logs).
        resp = self.session.post(f"{EUTILS}/{endpoint}", data=params, timeout=self.timeout)
        self._last_call = time.monotonic()
        resp.raise_for_status()
        return resp.json()

    def search_ids(self, query: str, max_ids: int = MAX_SEARCH_IDS) -> list[str]:
        data = self._eutils("esearch.fcgi", {"term": query, "retmax": max_ids})["esearchresult"]
        count = int(data.get("count", 0))
        if count > max_ids:
            logger.warning(
                "ClinVar returned %d records; only the first %d are ranked", count, max_ids
            )
        return list(data.get("idlist", []))

    def summaries(self, ids: Sequence[str]) -> list[dict[str, Any]]:
        docs: list[dict[str, Any]] = []
        for i in range(0, len(ids), ESUMMARY_BATCH):
            batch = ids[i : i + ESUMMARY_BATCH]
            result = self._eutils("esummary.fcgi", {"id": ",".join(batch)})["result"]
            docs.extend(result[uid] for uid in result.get("uids", []) if uid in result)
        return docs

    def anchor_base(self, chromosome: str, position: int) -> str:
        """GRCh38 reference base at a 1-based position (Ensembl)."""
        chrom = chromosome.removeprefix("chr")
        resp = self.session.get(
            f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{position}..{position}:1",
            headers={"Content-Type": "text/plain", "Accept": "text/plain"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        base = resp.text.strip().upper()
        if len(base) != 1:
            raise ValueError(f"unexpected reference sequence {base!r} at {chromosome}:{position}")
        return base

    def variants_for_gene(
        self,
        gene: str,
        *,
        consequence: str = "missense",
        significance: Sequence[str] = ("pathogenic", "likely_pathogenic"),
        max_variants: int | None = 25,
    ) -> list[ClinVarVariant]:
        """Known ClinVar variants for a gene, most trusted first."""
        query = build_query(gene, consequence, significance)
        logger.info("Searching ClinVar: %s", query)
        ids = self.search_ids(query)
        if not ids:
            return []
        parsed = [
            p
            for p in (parse_summary(d, gene) for d in self.summaries(ids))
            if p and matches_consequence(p, consequence)
        ]
        logger.info("%d ClinVar record(s), %d usable", len(ids), len(parsed))

        # Select before building variants so only the chosen indels cost an Ensembl call.
        selected: list[ClinVarVariant] = []
        for p in select_spread(parsed, max_variants):
            try:
                selected.append(self._to_variant(p))
            except (requests.RequestException, ValueError) as exc:
                logger.warning("Skipping %s: %s", p["accession"], exc)
        return selected

    def _to_variant(self, p: dict[str, Any]) -> ClinVarVariant:
        ref, alt, pos = p["ref"], p["alt"], p["position"]
        if not (ref and alt):  # pure indel: prepend the VCF anchor base
            pos -= 1
            anchor = self.anchor_base(p["chromosome"], pos)
            ref, alt = anchor + ref, anchor + alt
        return ClinVarVariant(
            variant=GenomicVariant(
                chromosome=p["chromosome"], position=pos, reference_bases=ref, alternate_bases=alt
            ),
            accession=p["accession"],
            title=p["title"],
            protein_change=p["protein_change"],
            classification=p["classification"],
            review_status=p["review_status"],
            consequences=p["consequences"],
        )


# --------------------------------------------------------------------------- VCF export


def write_vcf(variants: Iterable[ClinVarVariant], path: Path, source: str = "ClinVar") -> Path:
    """Write the selection as a VCF so the exact same set can be re-run with `vcf`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=geneme-structure {__version__} ({source})",
        "##reference=GRCh38",
        '##INFO=<ID=CLNVCV,Number=1,Type=String,Description="ClinVar VCV accession">',
        '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar germline classification">',
        '##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description="ClinVar review status">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]

    def info(value: str) -> str:
        return value.replace(" ", "_").replace(";", ",").replace("=", ":") or "."

    for cv in variants:
        v = cv.variant
        lines.append(
            "\t".join(
                [
                    v.chromosome,
                    str(v.position),
                    cv.label.replace(" ", "_"),
                    v.reference_bases,
                    v.alternate_bases,
                    ".",
                    ".",
                    f"CLNVCV={cv.accession};CLNSIG={info(cv.classification)};"
                    f"CLNREVSTAT={info(cv.review_status)}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")
    return path
