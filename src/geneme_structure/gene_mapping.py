"""Map genomic coordinates to genes and UniProt accessions.

* Variants -> Ensembl VEP REST (gene, consequence, protein position, UniProt)
* Regions  -> Ensembl overlap REST (protein-coding genes) -> UniProt REST search
"""

from __future__ import annotations

import logging
from typing import Any

import requests

from geneme_structure.http_session import build_session
from geneme_structure.models import GenomicInterval, GenomicVariant, VariantAnnotation

logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"
UNIPROT_REST = "https://rest.uniprot.org"
MAX_OVERLAP_BP = 5_000_000  # Ensembl overlap endpoint limit


class GeneMapper:
    def __init__(self, *, timeout: float = 30.0, session: requests.Session | None = None) -> None:
        self.timeout = timeout
        self.session = session or build_session()

    # ------------------------------------------------------------------ variants

    def annotate_variant(self, variant: GenomicVariant) -> VariantAnnotation:
        """Annotate a GRCh38 variant with Ensembl VEP (VCF-style input handles indels)."""
        vcf_line = (
            f"{variant.ensembl_chromosome} {variant.position} . "
            f"{variant.reference_bases} {variant.alternate_bases} . . ."
        )
        logger.info("Annotating %s with Ensembl VEP", variant)
        resp = self.session.post(
            f"{ENSEMBL_REST}/vep/homo_sapiens/region",
            json={"variants": [vcf_line], "uniprot": 1, "canonical": 1},
            headers={"Content-Type": "application/json"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        results = resp.json()
        if not results:
            return VariantAnnotation()
        return self._pick_consequence(results[0])

    @staticmethod
    def _pick_consequence(record: dict[str, Any]) -> VariantAnnotation:
        consequences: list[dict[str, Any]] = record.get("transcript_consequences") or []
        most_severe = record.get("most_severe_consequence")

        def rank(tc: dict[str, Any]) -> tuple[int, int, int]:
            # Prefer canonical transcripts that map to Swiss-Prot and have a protein position.
            return (
                int(bool(tc.get("canonical"))),
                int(bool(tc.get("swissprot"))),
                int(tc.get("protein_start") is not None),
            )

        if not consequences:
            return VariantAnnotation(most_severe_consequence=most_severe)

        best = max(consequences, key=rank)
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
        )

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
