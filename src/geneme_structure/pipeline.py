"""Genome -> regulatory prediction -> 3D structure pipeline, plus CLI.

Examples
--------
    # TP53 R175H hotspot (GRCh38), RNA-seq + DNase in transverse colon
    geneme-structure variant 'chr17:7675088:C>T' --condition UBERON:0001157

    # Compare tissues / cell types side by side
    geneme-structure variant 'chr17:7675088:C>T' \\
        --condition UBERON:0001157 --condition UBERON:0002107 --condition UBERON:0000178_CL:0001054

    # All protein-coding genes in a region
    geneme-structure region chr17:7668421-7687490 --condition UBERON:0001157

    # Batch: every variant in a VCF -> one summary table
    geneme-structure vcf examples/tp53_hotspots.vcf --max-variants 10

    # Start from a gene: known pathogenic ClinVar variants, most trusted first
    geneme-structure gene BRCA1 --max-variants 10 --consequence missense

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
import requests

from geneme_structure.alphafold_client import AlphaFoldClient
from geneme_structure.alphagenome_client import (
    SUPPORTED_SEQUENCE_LENGTHS,
    TRACK_OUTPUT_TYPES,
    AlphaGenomeClient,
    AlphaGenomeUnavailableError,
    ConditionLike,
    compare_conditions,
)
from geneme_structure.clinvar import (
    CONSEQUENCE_SETS,
    SIGNIFICANCE_FILTERS,
    ClinVarClient,
    ClinVarVariant,
    write_vcf,
)
from geneme_structure.config import MissingCredentialError, Settings
from geneme_structure.gene_mapping import GeneMapper
from geneme_structure.models import (
    GenomicInterval,
    GenomicVariant,
    PipelineResult,
    VariantAnnotation,
    VariantScores,
)
from geneme_structure.report import write_structure_html
from geneme_structure.vcf import VcfRecord, read_vcf

logger = logging.getLogger("geneme_structure")

DEFAULT_CONDITIONS: tuple[str, ...] = ("UBERON:0001157",)  # transverse colon
MAX_REGION_GENES = 10
MAX_VCF_VARIANTS = 25  # protects the AlphaGenome quota


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
        self.clinvar = ClinVarClient(
            timeout=settings.http_timeout,
            session=self.mapper.session,
            api_key=settings.ncbi_api_key.get_secret_value() if settings.ncbi_api_key else None,
        )
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
        conditions: Sequence[ConditionLike] = DEFAULT_CONDITIONS,
        output_types: Sequence[str] = (),
        uniprot_override: str | None = None,
        score: bool = True,
        annotation: VariantAnnotation | None = None,
    ) -> PipelineResult:
        """Annotate, score and fetch the structure for one variant.

        `score` runs AlphaGenome's recommended variant scorers (the headline
        numbers). `output_types` additionally requests raw REF/ALT tracks.
        Each condition is a separate (separately cached) AlphaGenome call.
        Pass `annotation` to reuse a batch VEP result.
        """
        warnings: list[str] = []

        if annotation is None:
            annotation = self.mapper.annotate_variant(variant)
        logger.info(
            "%s -> %s (%s), protein pos %s",
            variant,
            annotation.gene_symbol,
            annotation.most_severe_consequence,
            annotation.protein_position,
        )

        predictions, scores = [], []
        if self.alphagenome is not None:
            for condition in conditions:
                if score:
                    scores.append(self.alphagenome.score_variant(variant, condition))
                if output_types:
                    predictions.append(
                        self.alphagenome.predict_variant(variant, condition, output_types)
                    )

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
            alphagenome=predictions,
            variant_scores=scores,
            structures=structures,
            warnings=warnings,
        )

    def run_region(
        self,
        interval: GenomicInterval,
        conditions: Sequence[ConditionLike] = DEFAULT_CONDITIONS,
        output_types: Sequence[str] = ("RNA_SEQ",),
        max_genes: int = MAX_REGION_GENES,
    ) -> PipelineResult:
        warnings: list[str] = []
        genes = self.mapper.genes_in_region(interval)
        logger.info("%d protein-coding genes in %s", len(genes), interval)
        if len(genes) > max_genes:
            warnings.append(f"{len(genes)} genes found; fetching structures for first {max_genes}")

        predictions = []
        if self.alphagenome is not None:
            predictions = [
                self.alphagenome.predict_interval(interval, condition, output_types)
                for condition in conditions
            ]

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
            alphagenome=predictions,
            structures=structures,
            warnings=warnings,
        )

    def run_vcf(
        self,
        path: Path,
        conditions: Sequence[ConditionLike] = DEFAULT_CONDITIONS,
        max_variants: int | None = MAX_VCF_VARIANTS,
        score: bool = True,
    ) -> list[PipelineResult]:
        """Run every variant in a VCF. VEP annotation is batched; failures are per variant."""
        records = list(read_vcf(path, max_variants=max_variants))
        if not records:
            raise ValueError(f"No scorable variants in {path}")
        return self.run_records(records, conditions, score=score)

    def run_gene(
        self,
        gene: str,
        conditions: Sequence[ConditionLike] = DEFAULT_CONDITIONS,
        *,
        consequence: str = "missense",
        significance: Sequence[str] = ("pathogenic", "likely_pathogenic"),
        max_variants: int | None = MAX_VCF_VARIANTS,
        score: bool = True,
    ) -> tuple[list[ClinVarVariant], list[PipelineResult]]:
        """Known ClinVar variants for a gene -> the full pipeline for each.

        Returns the ClinVar selection too, so it can be saved as a VCF and re-run.
        """
        selection = self.clinvar.variants_for_gene(
            gene, consequence=consequence, significance=significance, max_variants=max_variants
        )
        if not selection:
            raise ValueError(
                f"No ClinVar variants for {gene!r} with consequence={consequence!r} and "
                f"significance={','.join(significance)}"
            )
        records = [
            VcfRecord(variant=cv.variant, vcf_id=cv.label, line_no=i)
            for i, cv in enumerate(selection, 1)
        ]
        return selection, self.run_records(records, conditions, score=score)

    def run_records(
        self,
        records: Sequence[VcfRecord],
        conditions: Sequence[ConditionLike] = DEFAULT_CONDITIONS,
        score: bool = True,
    ) -> list[PipelineResult]:
        """Run a batch of variants. VEP annotation is batched; failures are per variant."""
        annotations = self.mapper.annotate_variants([r.variant for r in records])

        results: list[PipelineResult] = []
        for i, (record, annotation) in enumerate(zip(records, annotations, strict=True), 1):
            logger.info("[%d/%d] %s", i, len(records), record.variant)
            try:
                result = self.run_variant(
                    record.variant, conditions, score=score, annotation=annotation
                )
            except (requests.RequestException, RuntimeError, ValueError) as exc:
                logger.error("%s failed: %s", record.variant, exc)
                result = PipelineResult(
                    query=str(record.variant), annotation=annotation, warnings=[f"failed: {exc}"]
                )
            if record.vcf_id:
                result = result.model_copy(update={"vcf_id": record.vcf_id})
            results.append(result)
        return results


# --------------------------------------------------------------------------- batch summary


def _top_effect(
    vs: VariantScores, output_types: set[str], scorer: str | None = None
) -> float | None:
    """Signed effect with the largest magnitude among the given output types (and scorer)."""
    hits = [
        s
        for s in vs.effects()
        if s.output_type in output_types and (scorer is None or s.scorer.startswith(scorer))
    ]
    return max(hits, key=lambda s: s.max_abs_score).top_score if hits else None


def batch_table(results: Sequence[PipelineResult]) -> pd.DataFrame:
    """One row per (variant, condition): annotation, ClinVar, structure and key effect scores.

    Score columns pick one representative scorer family each, so values are
    comparable down a column (same scorer), not across columns.
    """
    rows = []
    for r in results:
        a = r.annotation or VariantAnnotation()
        s = r.structures[0] if r.structures else None
        base = {
            "variant": r.query,
            "vcf_id": r.vcf_id,
            "rsid": a.rsid,
            "gene": a.gene_symbol,
            # All terms of the chosen transcript (VEP style), e.g. missense&splice_region
            "consequence": "&".join(a.consequence_terms) or a.most_severe_consequence,
            "protein_change": (
                f"{a.amino_acids}@{a.protein_position}" if a.protein_position else None
            ),
            "clinvar": ";".join(a.clinvar_significance) or None,
            "clinvar_somatic": a.clinvar_somatic,
            "gnomad_af": a.gnomad_af,
            "alphafold_entry": s.entry.entry_id if s else None,
            "residue_plddt": s.residue_plddt if s else None,
            "warnings": "; ".join(r.warnings) or None,
        }
        if not r.variant_scores:
            rows.append(base)
        for vs in r.variant_scores:
            rows.append(
                {
                    **base,
                    "condition": vs.label,
                    "rna_lfc": _top_effect(vs, {"RNA_SEQ"}, "GeneMaskLFC"),
                    "splicing": _top_effect(
                        vs, {"SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"}
                    ),
                    "accessibility": _top_effect(vs, {"DNASE", "ATAC"}),
                }
            )
    return pd.DataFrame(rows)


EFFECT_COLUMNS = {
    "rna_lfc": "expression effect (RNA-seq gene LFC)",
    "splicing": "splicing effect",
    "accessibility": "accessibility effect (DNase/ATAC)",
}


def batch_findings(table: pd.DataFrame, label_col: str = "vcf_id") -> list[str]:
    """Plain, data-driven sentences describing a `batch_table` (no biological claims).

    Used by the generic gene notebook so its conclusions follow the data for any gene.
    """
    if table.empty:
        return ["No variants were analysed."]
    labels = table[label_col].fillna(table["variant"]) if label_col in table else table["variant"]
    table = table.assign(_label=labels)
    per_variant = table.drop_duplicates("_label")
    out = [f"{len(per_variant)} variant(s) analysed."]

    plddt = per_variant["residue_plddt"].dropna() if "residue_plddt" in per_variant else []
    if len(plddt):
        bands = {
            "very high (>90)": int((plddt > 90).sum()),
            "confident (70-90)": int(((plddt > 70) & (plddt <= 90)).sum()),
            "low or very low (<=70)": int((plddt <= 70).sum()),
        }
        counts = ", ".join(f"{n} {band}" for band, n in bands.items() if n)
        out.append(
            f"AlphaFold confidence at the mutated residue ({len(plddt)} with a structure): "
            f"{counts}; median pLDDT {plddt.median():.1f}."
        )
    else:
        out.append("No AlphaFold structure covered these variants.")

    for col, name in EFFECT_COLUMNS.items():
        if col not in table or table[col].dropna().empty:
            continue
        # One condition at a time keeps scores comparable; use the first one.
        first = table[table["condition"] == table["condition"].dropna().iloc[0]]
        values = first.set_index("_label")[col].dropna()
        ranked = values.abs().sort_values(ascending=False)
        top = ranked.index[0]
        sentence = (
            f"Strongest {name}: {top} ({values[top]:+.3f}); "
            f"median |score| {ranked.median():.3f} across {len(ranked)} variant(s)"
        )
        if len(ranked) > 2 and ranked.iloc[0] > 2 * ranked.iloc[1]:
            second = ranked.index[1]
            sentence += f", more than twice the next ({second}, {values[second]:+.3f})"
        if col == "splicing" and "consequence" in first:
            terms = first.set_index("_label")["consequence"].fillna("")
            top3 = list(ranked.index[:3])
            near_sites = [lab for lab in top3 if "splice" in terms.get(lab, "")]
            if near_sites:
                sentence += (
                    f". Of the top {len(top3)}, {', '.join(near_sites)} "
                    f"{'is' if len(near_sites) == 1 else 'are'} annotated by VEP as "
                    "splice-region variant(s), an independent hint that the score is not noise"
                )
        out.append(sentence + ".")
    return out


# --------------------------------------------------------------------------- reporting


def report_stem(query: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in query).strip("_")


def write_report(result: PipelineResult, out_dir: Path) -> Path:
    """Write the JSON result plus CSV tables.

    * `<query>.json`: full result (per-track score records go to the CSV instead)
    * `<query>_scores.csv`: every AlphaGenome variant score (tidy, one row per track)
    * `<query>_conditions.csv`: effect scores side by side, when >1 condition
    * `<query>_tracks.csv`: raw-track summaries, when requested
    * `<query>_<AF entry>.html`: interactive 3D viewer for each fetched PDB structure
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = report_stem(result.query)
    json_path = out_dir / f"{stem}.json"
    json_path.write_text(
        result.model_dump_json(indent=2, exclude={"variant_scores": {"__all__": {"records"}}})
    )
    score_rows = [
        {"condition": vs.condition, **r.model_dump()}
        for vs in result.variant_scores
        for r in vs.records
    ]
    if score_rows:
        pd.DataFrame(score_rows).to_csv(out_dir / f"{stem}_scores.csv", index=False)
    if len(result.variant_scores) > 1:
        compare_conditions(result.variant_scores).to_csv(out_dir / f"{stem}_conditions.csv")
    track_rows = [
        {"condition": pred.condition, **t.model_dump()}
        for pred in result.alphagenome
        for t in pred.tracks
    ]
    if track_rows:
        pd.DataFrame(track_rows).to_csv(out_dir / f"{stem}_tracks.csv", index=False)
    for structure in result.structures:
        if structure.structure_path.endswith(".pdb"):
            html_path = out_dir / f"{stem}_{structure.entry.entry_id}.html"
            write_structure_html(result, structure, html_path)
            logger.info("3D viewer written to %s", html_path)
    return json_path


def print_summary(result: PipelineResult) -> None:
    print(f"\n=== {result.query} ===")
    if a := result.annotation:
        print(f"Gene:        {a.gene_symbol or '-'} ({a.gene_id or '-'})")
        print(f"Consequence: {a.most_severe_consequence or '-'}")
        if a.protein_position:
            print(f"Protein:     {a.amino_acids or '?'} at residue {a.protein_position}")
        if a.rsid or a.clinvar_significance or a.clinvar_somatic:
            sig = ", ".join(a.clinvar_significance) or "-"
            print(f"ClinVar:     {sig} (somatic: {a.clinvar_somatic or '-'}) {a.rsid or ''}")
        if a.gnomad_af is not None:
            print(f"gnomAD AF:   {a.gnomad_af:.3g}")
    if result.genes_in_region:
        print(f"Genes:       {', '.join(result.genes_in_region)}")
    fmt = {"float_format": lambda v: f"{v:+.4f}"}
    if len(result.variant_scores) == 1:
        vs = result.variant_scores[0]
        src = "cache" if vs.from_cache else "live"
        print(f"\nAlphaGenome variant scores ({src}, {vs.label} [{vs.condition}]):")
        df = pd.DataFrame([x.model_dump() for x in vs.effects()])
        if df.empty:
            print("  (no scores returned)")
        else:
            cols = ["output_type", "scorer", "n_scores", "top_score", "top_gene"]
            print(df[cols].fillna("-").to_string(index=False, **fmt))
    elif result.variant_scores:
        print("\nAlphaGenome effect scores by condition (compare within a row only):")
        for vs in result.variant_scores:
            src = "cache" if vs.from_cache else "live"
            print(f"  {vs.label:<40} {vs.condition} ({src})")
        print(compare_conditions(result.variant_scores).to_string(na_rep="-", **fmt))
    for pred in result.alphagenome:
        src = "cache" if pred.from_cache else "live"
        print(f"\nAlphaGenome raw tracks ({src}, {pred.condition}, window {pred.interval}):")
        df = pd.DataFrame([t.model_dump(exclude_none=True) for t in pred.tracks])
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
        "-c",
        "--condition",
        dest="conditions",
        action="append",
        metavar="CURIES",
        help="AlphaGenome condition: one ontology CURIE or several joined by '_' "
        "(max one per type), e.g. UBERON:0000178_CL:0001054. Repeat to compare "
        f"conditions (default: {DEFAULT_CONDITIONS[0]}, transverse colon)",
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

    pb = sub.add_parser("vcf", help="Analyse every variant in a GRCh38 VCF (.vcf or .vcf.gz)")
    pb.add_argument("vcf", type=Path)
    pb.add_argument(
        "--max-variants",
        type=int,
        default=MAX_VCF_VARIANTS,
        help="Stop after this many alleles; each costs AlphaGenome calls (default: %(default)s)",
    )
    pb.add_argument(
        "--reports", action="store_true", help="Also write per-variant JSON/CSV/HTML reports"
    )
    pb.add_argument(
        "--no-scores", action="store_true", help="Skip AlphaGenome recommended variant scorers"
    )
    _add_prediction_args(pb, default_outputs=[])

    pg = sub.add_parser(
        "gene", help="Analyse known pathogenic ClinVar variants of a gene, e.g. BRCA1"
    )
    pg.add_argument("gene", help="HGNC gene symbol")
    pg.add_argument(
        "--consequence",
        choices=sorted(CONSEQUENCE_SETS),
        default="missense",
        help="Which ClinVar variants to include (default: %(default)s)",
    )
    pg.add_argument(
        "--significance",
        default="pathogenic,likely_pathogenic",
        help=f"Comma-separated subset of {','.join(SIGNIFICANCE_FILTERS)} (default: %(default)s)",
    )
    pg.add_argument(
        "--max-variants",
        type=int,
        default=MAX_VCF_VARIANTS,
        help="Most trusted N variants (review status); each costs AlphaGenome calls "
        "(default: %(default)s)",
    )
    pg.add_argument(
        "--reports", action="store_true", help="Also write per-variant JSON/CSV/HTML reports"
    )
    pg.add_argument(
        "--no-scores", action="store_true", help="Skip AlphaGenome recommended variant scorers"
    )
    _add_prediction_args(pg, default_outputs=[])

    ps = sub.add_parser("structure", help="Fetch an AlphaFold DB model by UniProt accession")
    ps.add_argument("uniprot")
    ps.add_argument("--residue", type=int, help="Report pLDDT at this residue")
    ps.add_argument("--format", choices=("pdb", "cif"), default="pdb")
    return parser


def _run_vcf_command(pipeline: GenomeToStructurePipeline, args: argparse.Namespace) -> int:
    results = pipeline.run_vcf(
        args.vcf,
        args.conditions or DEFAULT_CONDITIONS,
        max_variants=args.max_variants,
        score=not args.no_scores,
    )
    return _write_batch(results, report_stem(args.vcf.name), args)


def _run_gene_command(pipeline: GenomeToStructurePipeline, args: argparse.Namespace) -> int:
    significance = [s.strip() for s in args.significance.split(",") if s.strip()]
    selection, results = pipeline.run_gene(
        args.gene,
        args.conditions or DEFAULT_CONDITIONS,
        consequence=args.consequence,
        significance=significance,
        max_variants=args.max_variants,
        score=not args.no_scores,
    )
    stem = f"{report_stem(args.gene)}_clinvar_{args.consequence}"
    vcf_path = write_vcf(selection, args.out_dir / f"{stem}.vcf")
    logger.info("ClinVar selection saved to %s (re-run it with the `vcf` command)", vcf_path)
    return _write_batch(results, stem, args)


def _write_batch(results: list[PipelineResult], stem: str, args: argparse.Namespace) -> int:
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.reports:
        for result in results:
            write_report(result, out_dir / stem)
    table = batch_table(results)
    csv_path = out_dir / f"{stem}_summary.csv"
    table.to_csv(csv_path, index=False)

    cols = [
        c
        for c in (
            "vcf_id",
            "variant",
            "gene",
            "protein_change",
            "clinvar",
            "residue_plddt",
            "condition",
            "rna_lfc",
            "splicing",
            "accessibility",
        )
        if c in table
    ]
    with pd.option_context("display.width", 200, "display.max_colwidth", 28):
        print(table[cols].to_string(index=False, na_rep="-", float_format=lambda v: f"{v:.4g}"))
    failed = sum(1 for r in results if any(w.startswith("failed") for w in r.warnings))
    logger.info("%d variant(s), %d failed. Summary written to %s", len(results), failed, csv_path)
    return 1 if failed == len(results) else 0


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
            if args.command == "vcf":
                return _run_vcf_command(pipeline, args)
            if args.command == "gene":
                return _run_gene_command(pipeline, args)
            if args.command == "variant":
                result = pipeline.run_variant(
                    args.variant,
                    args.conditions or DEFAULT_CONDITIONS,
                    args.outputs,
                    uniprot_override=args.uniprot,
                    score=not args.no_scores,
                )
            else:
                result = pipeline.run_region(
                    args.region,
                    args.conditions or DEFAULT_CONDITIONS,
                    args.outputs,
                    max_genes=args.max_genes,
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
