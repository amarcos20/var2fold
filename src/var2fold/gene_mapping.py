"""Map genomic coordinates to genes and UniProt accessions.

* Variants -> Ensembl VEP REST (gene, consequence, protein position, UniProt,
  plus ClinVar significance, rsID and gnomAD frequency from co-located variants)
* Regions  -> Ensembl overlap REST (protein-coding genes) -> UniProt REST search
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

import requests

from var2fold.http_session import build_session
from var2fold.models import (
    GenomicInterval,
    GenomicVariant,
    ProteinFeature,
    VariantAnnotation,
)

logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"
UNIPROT_REST = "https://rest.uniprot.org"
MAX_OVERLAP_BP = 5_000_000  # Ensembl overlap endpoint limit
VEP_BATCH_SIZE = 200  # Ensembl VEP POST limit
VEP_OPTIONS = {"uniprot": 1, "canonical": 1}
# UniProt feature types that describe domain architecture (Region is filtered further).
DOMAIN_FEATURE_TYPES = ("Domain", "DNA binding", "Zinc finger", "Region")
_REGION_NOISE = ("interaction", "required for", "disordered", "necessary for")


class GeneMapper:
    def __init__(self, *, timeout: float = 30.0, session: requests.Session | None = None) -> None:
        self.timeout = timeout
        self.session = session or build_session()

    # ------------------------------------------------------------------ variants

    def annotate_variant(self, variant: GenomicVariant) -> VariantAnnotation:
        """Annotate one GRCh38 variant with Ensembl VEP."""
        return self.annotate_variants([variant])[0]

    def annotate_variants(self, variants: Sequence[GenomicVariant]) -> list[VariantAnnotation]:
        """Annotate GRCh38 variants with Ensembl VEP, in batches of up to 200.

        VCF-style input handles indels correctly. Returns one annotation per
        input variant, in order (an empty annotation if VEP returns nothing).
        """
        out: list[VariantAnnotation] = []
        for i in range(0, len(variants), VEP_BATCH_SIZE):
            batch = list(variants[i : i + VEP_BATCH_SIZE])
            lines = [vep_input_line(v) for v in batch]
            logger.info("Annotating %d variant(s) with Ensembl VEP", len(batch))
            resp = self.session.post(
                f"{ENSEMBL_REST}/vep/homo_sapiens/region",
                json={"variants": lines, **VEP_OPTIONS},
                headers={"Content-Type": "application/json"},
                timeout=self.timeout * (1 + len(batch) // 50),
            )
            resp.raise_for_status()
            by_input = {rec.get("input"): rec for rec in resp.json()}
            for variant, line in zip(batch, lines, strict=True):
                record = by_input.get(line)
                out.append(
                    parse_vep_record(record, variant.alternate_bases)
                    if record
                    else VariantAnnotation()
                )
        return out

    # ------------------------------------------------------------------ regions

    def genes_in_region(
        self, interval: GenomicInterval, protein_coding_only: bool = True
    ) -> list[str]:
        if interval.width > MAX_OVERLAP_BP:
            raise ValueError(f"Region wider than {MAX_OVERLAP_BP:,} bp is not supported by Ensembl")
        chrom = interval.chromosome.removeprefix("chr").replace("M", "MT")
        region = f"{chrom}:{interval.start + 1}-{interval.end}"
        logger.info("Listing genes overlapping %s", interval)
        resp = self.session.get(
            f"{ENSEMBL_REST}/overlap/region/homo_sapiens/{region}",
            params={"feature": "gene"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        genes = [
            g.get("external_name") or g["id"]
            for g in resp.json()
            if not protein_coding_only or g.get("biotype") == "protein_coding"
        ]
        return list(dict.fromkeys(genes))

    def uniprot_for_gene(self, gene_symbol: str, taxon_id: int = 9606) -> str | None:
        """Reviewed (Swiss-Prot) accession for a gene symbol, or None."""
        resp = self.session.get(
            f"{UNIPROT_REST}/uniprotkb/search",
            params={
                "query": f"gene_exact:{gene_symbol} AND organism_id:{taxon_id} AND reviewed:true",
                "fields": "accession,gene_primary",
                "format": "json",
                "size": 1,
            },
            timeout=self.timeout,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])
        return results[0]["primaryAccession"] if results else None

    def protein_features(self, accession: str) -> list[ProteinFeature]:
        """Domain architecture of a UniProt entry (domains, DNA-binding, zinc fingers...)."""
        resp = self.session.get(
            f"{UNIPROT_REST}/uniprotkb/{accession}.json",
            params={"fields": "ft_domain,ft_region,ft_dna_bind,ft_zn_fing"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        return parse_uniprot_features(resp.json())


def parse_uniprot_features(entry: dict[str, Any]) -> list[ProteinFeature]:
    """Keep architectural features; drop interaction/disorder 'Region' annotations.

    Disorder is left out on purpose: the pLDDT profile already shows it.
    """
    out: list[ProteinFeature] = []
    for f in entry.get("features") or []:
        kind = f.get("type")
        if kind not in DOMAIN_FEATURE_TYPES:
            continue
        name = (f.get("description") or "").strip() or kind
        if kind == "Region" and name.lower().startswith(_REGION_NOISE):
            continue
        loc = f.get("location") or {}
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if isinstance(start, int) and isinstance(end, int) and end >= start:
            out.append(ProteinFeature(start=start, end=end, name=name, kind=kind))
    return sorted(out, key=lambda x: (x.start, x.end))


# --------------------------------------------------------------------------- VEP parsing


def vep_input_line(variant: GenomicVariant) -> str:
    """VEP's whitespace-separated VCF input; also used to match results back to inputs."""
    return (
        f"{variant.ensembl_chromosome} {variant.position} . "
        f"{variant.reference_bases} {variant.alternate_bases} . . ."
    )


def _rank_transcript(tc: dict[str, Any]) -> tuple[int, int, int]:
    # Prefer canonical transcripts that map to Swiss-Prot and have a protein position.
    return (
        int(bool(tc.get("canonical"))),
        int(bool(tc.get("swissprot"))),
        int(tc.get("protein_start") is not None),
    )


def _vep_alt(record: dict[str, Any], alt: str) -> str:
    """ALT allele as VEP reports it (indels lose their anchor base; deletions become '-')."""
    allele = record.get("allele_string", "")
    return allele.split("/", 1)[1] if "/" in allele else alt


def parse_clinvar(colocated: list[dict[str, Any]], alt: str) -> dict[str, Any]:
    """rsID, allele-specific ClinVar significance, somatic class, VCV ids and gnomAD AF."""
    rsid: str | None = None
    somatic: str | None = None
    significance: list[str] = []
    vcv_ids: set[str] = set()
    gnomad: list[float] = []
    # ClinVar writes a deleted allele as '' where VEP uses '-'
    alt_forms = {alt, ""} if alt == "-" else {alt}

    for cv in colocated:
        vid = str(cv.get("id", ""))
        if rsid is None and vid.startswith("rs"):
            rsid = vid
        if somatic is None and cv.get("clinvar_somatic_classification"):
            somatic = cv["clinvar_somatic_classification"]
        # 'A:likely_pathogenic;T:pathogenic/likely_pathogenic' -> our ALT only
        for part in (cv.get("clin_sig_allele") or "").split(";"):
            allele, _, sig = part.partition(":")
            if allele in alt_forms and sig:
                significance.extend(sig.split("/"))
        synonyms = cv.get("var_synonyms")
        if isinstance(synonyms, dict):
            vcv_ids.update(str(x) for x in synonyms.get("ClinVar", []) if str(x).startswith("VCV"))
        freqs = (cv.get("frequencies") or {}).get(alt) or {}
        gnomad.extend(
            float(freqs[k]) for k in ("gnomadg", "gnomade") if isinstance(freqs.get(k), int | float)
        )

    return {
        "rsid": rsid,
        "clinvar_significance": list(dict.fromkeys(significance)),
        "clinvar_somatic": somatic,
        "clinvar_ids": sorted(vcv_ids),
        "gnomad_af": max(gnomad) if gnomad else None,
    }


def parse_vep_record(record: dict[str, Any], alt: str) -> VariantAnnotation:
    """Turn one VEP JSON record into a `VariantAnnotation`."""
    clinvar = parse_clinvar(record.get("colocated_variants") or [], _vep_alt(record, alt))
    consequences: list[dict[str, Any]] = record.get("transcript_consequences") or []
    most_severe = record.get("most_severe_consequence")
    if not consequences:
        return VariantAnnotation(most_severe_consequence=most_severe, **clinvar)

    best = max(consequences, key=_rank_transcript)
    accessions = [acc.split(".")[0] for acc in best.get("swissprot") or []]
    return VariantAnnotation(
        gene_symbol=best.get("gene_symbol"),
        gene_id=best.get("gene_id"),
        transcript_id=best.get("transcript_id"),
        most_severe_consequence=most_severe,
        consequence_terms=best.get("consequence_terms") or [],
        protein_position=best.get("protein_start"),
        amino_acids=best.get("amino_acids"),
        uniprot_accessions=list(dict.fromkeys(accessions)),
        **clinvar,
    )
