"""Genome -> regulatory prediction -> 3D structure pipeline, plus CLI.

Examples
--------
    # TP53 R175H hotspot (GRCh38), RNA-seq + DNase in transverse colon
    geneme-structure variant 'chr17:7675088:C>T' --ontology UBERON:0001157 --outputs RNA_SEQ DNASE

    # All protein-coding genes in a region
    geneme-structure region chr17:7668421-7687490 --ontology UBERON:0001157

    # Structure only (no AlphaGenome key needed)
    geneme-structure structure P04637 --residue 175

`python -m geneme_structure.pipeline ...` works too.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from geneme_structure.alphafold_client import AlphaFoldClient
from geneme_structure.alphagenome_client import (
    SUPPORTED_SEQUENCE_LENGTHS,
    TRACK_OUTPUT_TYPES,
    AlphaGenomeClient,
    AlphaGenomeUnavailableError,
)
from geneme_structure.config import MissingCredentialError, Settings
from geneme_structure.gene_mapping import GeneMapper
from geneme_structure.models import GenomicInterval, GenomicVariant, PipelineResult

logger = logging.getLogger("geneme_structure")

DEFAULT_ONTOLOGY = ("UBERON:0001157",)  # transverse colon
MAX_REGION_GENES = 10


class GenomeToStructurePipeline:
    """Orchestrates annotation (Ensembl), AlphaGenome and AlphaFold DB."""

    def __init__(
        self,
        settings: Settings,
        *,
        use_alphagenome: bool = True,
        use_cache: bool = True,
        sequence_length: int = 1_048_576,
    ) -> None:
        self.settings = settings
        self.mapper = GeneMapper(timeout=settings.http_timeout)
        self.alphafold = AlphaFoldClient(settings.data_dir, timeout=settings.http_timeout)
        self.alphagenome: AlphaGenomeClient | None = None
        if use_alphagenome:
            self.alphagenome = AlphaGenomeClient(
                settings.require_alphagenome_key(),
                settings.cache_dir,
                sequence_length=sequence_length,
                use_cache=use_cache,
            )

    def run_variant(
        self,
        variant: GenomicVariant,
        ontology_terms: Sequence[str] = DEFAULT_ONTOLOGY,
        output_types: Sequence[str] = (),
        uniprot_override: str | None = None,
        score: bool = True,
    ) -> PipelineResult:
        """Annotate, score and fetch the structure for one variant.

        `score` runs AlphaGenome's recommended variant scorers (the headline
        numbers). `output_types` additionally requests raw REF/ALT tracks.
        """
        warnings: list[str] = []

        annotation = self.mapper.annotate_variant(variant)
        logger.info(
            "%s -> %s (%s), protein pos %s",
            variant,
            annotation.gene_symbol,
            annotation.most_severe_consequence,
            annotation.protein_position,
        )

        prediction = scores = None
        if self.alphagenome is not None:
            if score:
                scores = self.alphagenome.score_variant(variant, ontology_terms)
            if output_types:
                prediction = self.alphagenome.predict_variant(variant, ontology_terms, output_types)

        accession = uniprot_override or annotation.primary_uniprot
        if accession is None and annotation.gene_symbol:
            accession = self.mapper.uniprot_for_gene(annotation.gene_symbol)

        structures = []
        if accession:
            residue = annotation.protein_position if not uniprot_override else None
            if s := self.alphafold.fetch_structure(accession, residue=residue):
                structures.append(s)
            else:
                warnings.append(f"No AlphaFold DB model for {accession}")
        else:
            warnings.append("Variant does not map to a protein-coding gene; no structure fetched")

        return PipelineResult(
            query=str(variant),
            annotation=annotation,
            alphagenome=prediction,
            variant_scores=scores,
            structures=structures,
            warnings=warnings,
        )

    def run_region(
        self,
        interval: GenomicInterval,
        ontology_terms: Sequence[str] = DEFAULT_ONTOLOGY,
        output_types: Sequence[str] = ("RNA_SEQ",),
        max_genes: int = MAX_REGION_GENES,
    ) -> PipelineResult:
        warnings: list[str] = []
        genes = self.mapper.genes_in_region(interval)
        logger.info("%d protein-coding genes in %s", len(genes), interval)
        if len(genes) > max_genes:
            warnings.append(f"{len(genes)} genes found; fetching structures for first {max_genes}")

        prediction = None
        if self.alphagenome is not None:
            prediction = self.alphagenome.predict_interval(interval, ontology_terms, output_types)

        structures = []
        for gene in genes[:max_genes]:
            accession = self.mapper.uniprot_for_gene(gene)
            if accession is None:
                warnings.append(f"No reviewed UniProt entry for {gene}")
                continue
            if s := self.alphafold.fetch_structure(accession):
                structures.append(s)
            else:
                warnings.append(f"No AlphaFold DB model for {gene} ({accession})")

        return PipelineResult(
            query=str(interval),
            genes_in_region=genes,
            alphagenome=prediction,
            structures=structures,
            warnings=warnings,
        )


# --------------------------------------------------------------------------- reporting


def report_stem(query: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in query).strip("_")


def write_report(result: PipelineResult, out_dir: Path) -> Path:
    """Write the JSON result plus CSV tables.

    * `<query>.json`: full result (per-track score records go to the CSV instead)
    * `<query>_scores.csv`: every AlphaGenome variant score (tidy, one row per track)
    * `<query>_tracks.csv`: raw-track summaries, when requested
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = report_stem(result.query)
    json_path = out_dir / f"{stem}.json"
    json_path.write_text(result.model_dump_json(indent=2, exclude={"variant_scores": {"records"}}))
    if result.variant_scores and result.variant_scores.records:
        df = pd.DataFrame([r.model_dump() for r in result.variant_scores.records])
        df.to_csv(out_dir / f"{stem}_scores.csv", index=False)
    if result.alphagenome and result.alphagenome.tracks:
        df = pd.DataFrame([t.model_dump() for t in result.alphagenome.tracks])
        df.to_csv(out_dir / f"{stem}_tracks.csv", index=False)
    return json_path


def print_summary(result: PipelineResult) -> None:
    print(f"\n=== {result.query} ===")
    if a := result.annotation:
        print(f"Gene:        {a.gene_symbol or '-'} ({a.gene_id or '-'})")
        print(f"Consequence: {a.most_severe_consequence or '-'}")
        if a.protein_position:
            print(f"Protein:     {a.amino_acids or '?'} at residue {a.protein_position}")
    if result.genes_in_region:
        print(f"Genes:       {', '.join(result.genes_in_region)}")
    if vs := result.variant_scores:
        src = "cache" if vs.from_cache else "live"
        print(f"\nAlphaGenome variant scores ({src}, {vs.condition}, window {vs.interval}):")
        df = pd.DataFrame([x.model_dump() for x in vs.effects()])
        if df.empty:
            print("  (no scores returned)")
        else:
            cols = ["output_type", "scorer", "n_scores", "top_score", "top_gene"]
            print(df[cols].fillna("-").to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    if p := result.alphagenome:
        src = "cache" if p.from_cache else "live"
        print(f"\nAlphaGenome ({src}, window {p.interval}, {', '.join(p.ontology_terms)}):")
        df = pd.DataFrame([t.model_dump(exclude_none=True) for t in p.tracks])
        print(df.to_string(index=False) if not df.empty else "  (no tracks returned)")
    for s in result.structures:
        e = s.entry
        print(f"\nAlphaFold:   {e.entry_id} - {e.description or e.gene or ''}")
        print(f"  file:      {s.structure_path}")
        print(f"  mean pLDDT {s.mean_plddt}")
        if s.residue_plddt is not None:
            print(f"  residue {s.residue}: pLDDT {s.residue_plddt:.1f} ({s.residue_confidence})")
    for w in result.warnings:
        print(f"! {w}")


# --------------------------------------------------------------------------- CLI


def _add_prediction_args(p: argparse.ArgumentParser, default_outputs: list[str]) -> None:
    p.add_argument(
        "--ontology",
        nargs="+",
        default=list(DEFAULT_ONTOLOGY),
        help="Ontology CURIEs combined into one AlphaGenome condition, max one per type, "
        "e.g. UBERON:0000178 CL:0001054 (default: %(default)s)",
    )
    p.add_argument(
        "--outputs",
        nargs="+",
        default=default_outputs,
        choices=TRACK_OUTPUT_TYPES,
        metavar="TYPE",
        help=f"Raw AlphaGenome tracks to summarise: {', '.join(TRACK_OUTPUT_TYPES)} "
        "(default: %(default)s)",
    )
    p.add_argument(
        "--sequence-length",
        type=int,
        default=1_048_576,
        choices=SUPPORTED_SEQUENCE_LENGTHS,
        help="Context window in bp (default: %(default)s)",
    )
    p.add_argument("--no-alphagenome", action="store_true", help="Skip AlphaGenome predictions")
    p.add_argument("--no-cache", action="store_true", help="Bypass the AlphaGenome response cache")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="geneme-structure",
        description="Link variants/regions to AlphaGenome predictions and AlphaFold DB structures.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    parser.add_argument("-o", "--out-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    pv = sub.add_parser("variant", help="Analyse a single GRCh38 variant, e.g. chr17:7675088:C>T")
    pv.add_argument("variant", type=GenomicVariant.parse)
    pv.add_argument("--uniprot", help="Override the UniProt accession used for the structure")
    pv.add_argument(
        "--no-scores", action="store_true", help="Skip AlphaGenome recommended variant scorers"
    )
    _add_prediction_args(pv, default_outputs=[])

    pr = sub.add_parser("region", help="Analyse a GRCh38 region, e.g. chr17:7668421-7687490")
    pr.add_argument("region", type=GenomicInterval.parse)
    pr.add_argument("--max-genes", type=int, default=MAX_REGION_GENES)
    _add_prediction_args(pr, default_outputs=["RNA_SEQ"])

    ps = sub.add_parser("structure", help="Fetch an AlphaFold DB model by UniProt accession")
    ps.add_argument("uniprot")
    ps.add_argument("--residue", type=int, help="Report pLDDT at this residue")
    ps.add_argument("--format", choices=("pdb", "cif"), default="pdb")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    # Silence noisy third-party loggers unless debugging.
    for noisy in ("urllib3", "absl", "grpc"):
        logging.getLogger(noisy).setLevel(logging.DEBUG if args.verbose else logging.WARNING)

    settings = Settings.from_env()

    try:
        if args.command == "structure":
            client = AlphaFoldClient(settings.data_dir, timeout=settings.http_timeout)
            s = client.fetch_structure(args.uniprot, residue=args.residue, fmt=args.format)
            if s is None:
                logger.error("No AlphaFold DB model for %s", args.uniprot)
                return 1
            result = PipelineResult(query=args.uniprot, structures=[s])
        else:
            pipeline = GenomeToStructurePipeline(
                settings,
                use_alphagenome=not args.no_alphagenome,
                use_cache=not args.no_cache,
                sequence_length=args.sequence_length,
            )
            if args.command == "variant":
                result = pipeline.run_variant(
                    args.variant,
                    args.ontology,
                    args.outputs,
                    uniprot_override=args.uniprot,
                    score=not args.no_scores,
                )
            else:
                result = pipeline.run_region(
                    args.region, args.ontology, args.outputs, max_genes=args.max_genes
                )
    except (MissingCredentialError, AlphaGenomeUnavailableError) as exc:
        logger.error("%s (or pass --no-alphagenome)", exc)
        return 2
    except ValueError as exc:
        logger.error("%s", exc)
        return 2

    print_summary(result)
    path = write_report(result, args.out_dir)
    logger.info("Report written to %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
