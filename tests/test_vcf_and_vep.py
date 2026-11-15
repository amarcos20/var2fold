import gzip
from pathlib import Path

import pytest

from geneme_structure.gene_mapping import parse_clinvar, parse_vep_record, vep_input_line
from geneme_structure.models import (
    GenomicInterval,
    GenomicVariant,
    PipelineResult,
    ScoreRecord,
    VariantAnnotation,
    VariantScores,
)
from geneme_structure.pipeline import batch_table
from geneme_structure.vcf import read_vcf

VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr17\t7675088\tR175H\tC\tT\t.\t.\t.
17\t100\t.\tA\tG,T\t.\t.\t.
chr1\t200\tsv\tN\t<DEL>\t.\t.\t.
chr1\t300\tstar\tA\t*\t.\t.\t.
chrUn_KI270302v1\t5\tbadcontig\tA\tG\t.\t.\t.
chr2\t400\tbadbase\tA\tR\t.\t.\t.
"""


@pytest.mark.parametrize("compress", [False, True])
def test_read_vcf_splits_and_skips(tmp_path: Path, compress: bool):
    path = tmp_path / ("x.vcf.gz" if compress else "x.vcf")
    if compress:
        with gzip.open(path, "wt") as fh:
            fh.write(VCF)
    else:
        path.write_text(VCF)
    records = list(read_vcf(path))
    assert [str(r.variant) for r in records] == [
        "chr17:7675088:C>T",
        "chr17:100:A>G",
        "chr17:100:A>T",
    ]
    assert records[0].vcf_id == "R175H" and records[1].vcf_id is None
    assert len(list(read_vcf(path, max_variants=2))) == 2


def test_read_vcf_rejects_malformed(tmp_path: Path):
    path = tmp_path / "bad.vcf"
    path.write_text("chr1 100 . A G\n")  # spaces, not tabs
    with pytest.raises(ValueError):
        list(read_vcf(path))


def test_vep_input_line_uses_ensembl_contig():
    assert vep_input_line(GenomicVariant.parse("chrM:10:A>G")) == "MT 10 . A G . . ."


def test_parse_clinvar_is_allele_specific():
    colocated = [
        {"id": "COSV1", "clinvar_somatic_classification": "Neoplasm(oncogenicity:Oncogenic)"},
        {
            "id": "rs28934578",
            "clin_sig_allele": "A:likely_pathogenic;T:pathogenic/likely_pathogenic",
            "var_synonyms": {"ClinVar": ["RCV1", "VCV000012374"], "OMIM": [1]},
            "frequencies": {"T": {"gnomadg": 1e-5, "gnomade": 6e-6, "gnomade_nfe": 0.1}},
        },
    ]
    info = parse_clinvar(colocated, "T")
    assert info == {
        "rsid": "rs28934578",
        "clinvar_significance": ["pathogenic", "likely_pathogenic"],
        "clinvar_somatic": "Neoplasm(oncogenicity:Oncogenic)",
        "clinvar_ids": ["VCV000012374"],
        "gnomad_af": 1e-5,  # overall genome/exome only, not sub-populations
    }
    assert parse_clinvar(colocated, "A")["clinvar_significance"] == ["likely_pathogenic"]


def test_parse_clinvar_deletion_uses_empty_allele():
    colocated = [{"id": "rs113993960", "clin_sig_allele": ":pathogenic;:drug_response"}]
    assert parse_clinvar(colocated, "-")["clinvar_significance"] == ["pathogenic", "drug_response"]


def test_parse_vep_record_prefers_canonical_swissprot():
    record = {
        "allele_string": "C/T",
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [
            {"gene_symbol": "TP53", "transcript_id": "ENST_other", "protein_start": 43},
            {
                "gene_symbol": "TP53",
                "transcript_id": "ENST00000269305",
                "canonical": 1,
                "swissprot": ["P04637.250"],
                "protein_start": 175,
                "amino_acids": "R/H",
            },
        ],
    }
    a = parse_vep_record(record, "T")
    assert (a.transcript_id, a.protein_position, a.primary_uniprot) == (
        "ENST00000269305",
        175,
        "P04637",
    )


def test_batch_table_one_row_per_variant_and_condition():
    def scores(label: str, lfc: float) -> VariantScores:
        from geneme_structure.alphagenome_client import summarize_scores

        recs = [
            ScoreRecord(
                output_type="RNA_SEQ",
                variant_scorer="GeneMaskLFCScorer(DIFF_MEAN)",
                biosample_name=label,
                raw_score=lfc,
            ),
            ScoreRecord(
                output_type="SPLICE_SITES",
                variant_scorer="GeneMaskSplicingScorer()",
                biosample_name=label,
                raw_score=0.3,
            ),
        ]
        return VariantScores(
            query="v",
            interval=GenomicInterval(chromosome="chr1", start=0, end=16_384),
            condition=label,
            summary=summarize_scores(recs),
            records=recs,
        )

    annotated = PipelineResult(
        query="chr17:7675088:C>T",
        vcf_id="R175H",
        annotation=VariantAnnotation(
            gene_symbol="TP53",
            amino_acids="R/H",
            protein_position=175,
            clinvar_significance=["pathogenic"],
        ),
        variant_scores=[scores("colon", -0.1), scores("liver", 0.2)],
    )
    bare = PipelineResult(query="chr1:5:A>G", warnings=["failed: boom"])
    table = batch_table([annotated, bare])
    assert len(table) == 3
    first = table.iloc[0]
    assert (first.vcf_id, first.protein_change, first.clinvar) == (
        "R175H",
        "R/H@175",
        "pathogenic (submissions)",
    )
    assert list(table.condition[:2]) == ["colon", "liver"]
    assert list(table.rna_lfc[:2]) == [-0.1, 0.2] and table.splicing[0] == 0.3
    assert table.iloc[2].warnings == "failed: boom"
