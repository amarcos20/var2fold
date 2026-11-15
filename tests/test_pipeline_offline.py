"""End-to-end pipeline and CLI tests with every network dependency faked."""

import json
from pathlib import Path

import pandas as pd
import pytest
from fakes import AF_ENTRY, TOY_PDB, FakeResponse, FakeSession, clinvar_doc, eutils, vep_echo

from geneme_structure import alphafold_client, clinvar, gene_mapping
from geneme_structure.alphagenome_client import AlphaGenomeClient, summarize_scores
from geneme_structure.config import Settings
from geneme_structure.models import GenomicVariant, ScoreRecord, VariantScores
from geneme_structure.pipeline import GenomeToStructurePipeline, main, write_report

VEP_TP53 = {
    "allele_string": "C/T",
    "most_severe_consequence": "missense_variant",
    "transcript_consequences": [
        {
            "gene_symbol": "TP53",
            "gene_id": "ENSG00000141510",
            "canonical": 1,
            "swissprot": ["P04637"],
            "protein_start": 2,
            "amino_acids": "R/H",
        }
    ],
    "colocated_variants": [{"id": "rs1", "clin_sig_allele": "T:pathogenic"}],
}


def _session() -> FakeSession:
    s = FakeSession()
    s.add("/vep/", vep_echo(lambda line: VEP_TP53))
    s.add("/api/prediction/", FakeResponse(payload=[AF_ENTRY]))
    s.add("/files/", FakeResponse(content=TOY_PDB.encode()))
    s.add("/info/software", FakeResponse(payload={"release": 116}))
    return s


class FakeAlphaGenome:
    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []

    def score_variant(self, variant: GenomicVariant, condition: str) -> VariantScores:
        self.calls.append((str(variant), condition))
        recs = [
            ScoreRecord(
                output_type="RNA_SEQ",
                variant_scorer="GeneMaskLFCScorer(DIFF_MEAN)",
                gene_name="TP53",
                biosample_name=f"tissue {condition}",
                raw_score=-0.25,
            )
        ]
        window = variant.reference_interval.resized(16_384)
        return VariantScores(
            query=str(variant),
            interval=window,
            condition=condition,
            summary=summarize_scores(recs),
            records=recs,
        )


@pytest.fixture
def settings(tmp_path: Path) -> Settings:
    return Settings(cache_dir=tmp_path / "cache", data_dir=tmp_path / "data")


def test_run_variant_end_to_end(settings, tmp_path):
    pipeline = GenomeToStructurePipeline(settings, use_alphagenome=False)
    session = _session()
    pipeline.mapper = gene_mapping.GeneMapper(session=session)
    pipeline.alphafold = alphafold_client.AlphaFoldClient(settings.data_dir, session=session)
    pipeline.alphagenome = FakeAlphaGenome()  # type: ignore[assignment]

    result = pipeline.run_variant(
        GenomicVariant.parse("chr17:7675088:C>T"), conditions=["UBERON:1", "UBERON:2"]
    )
    assert result.annotation.clinvar_significance == ["pathogenic"]
    assert [vs.condition for vs in result.variant_scores] == ["UBERON:1", "UBERON:2"]
    assert result.structures[0].residue_plddt == 72.5
    assert result.warnings == []

    json_path = write_report(result, tmp_path / "out")
    produced = sorted(p.name for p in (tmp_path / "out").iterdir())
    assert produced == [
        "chr17_7675088_C_T.json",
        "chr17_7675088_C_T_AF-P04637-F1.html",
        "chr17_7675088_C_T_conditions.csv",
        "chr17_7675088_C_T_scores.csv",
    ]
    saved = json.loads(json_path.read_text())
    assert "records" not in saved["variant_scores"][0]  # per-track rows live in the CSV
    scores = pd.read_csv(tmp_path / "out" / "chr17_7675088_C_T_scores.csv")
    assert set(scores.condition) == {"UBERON:1", "UBERON:2"}


def test_non_coding_variant_warns_instead_of_failing(settings):
    pipeline = GenomeToStructurePipeline(settings, use_alphagenome=False)
    session = FakeSession().add(
        "/vep/", vep_echo(lambda line: {"most_severe_consequence": "intergenic_variant"})
    )
    pipeline.mapper = gene_mapping.GeneMapper(session=session)
    result = pipeline.run_variant(GenomicVariant.parse("chr1:1000:A>G"))
    assert result.structures == []
    assert "no structure fetched" in result.warnings[0]


def test_cli_vcf_without_alphagenome(tmp_path, monkeypatch, capsys):
    session = _session()
    monkeypatch.setattr(gene_mapping, "build_session", lambda: session)
    monkeypatch.setattr(alphafold_client, "build_session", lambda: session)
    monkeypatch.chdir(tmp_path)
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\nchr17\t7675088\tR175H\tC\tT\nchr17\t7674220\tR248Q\tC\tT\n"
    )

    code = main(["-o", str(tmp_path / "res"), "vcf", str(vcf), "--no-alphagenome"])
    assert code == 0
    table = pd.read_csv(tmp_path / "res" / "in_vcf_summary.csv")
    assert list(table.vcf_id) == ["R175H", "R248Q"]
    assert list(table.residue_plddt) == [72.5, 72.5]
    vep_posts = [c for c in session.calls if "/vep/" in c[1]]
    assert len(vep_posts) == 1  # both variants annotated in one batch
    assert "R/H@2" in capsys.readouterr().out


def test_cli_reports_missing_key(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # no .env here
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    assert main(["variant", "chr17:7675088:C>T"]) == 2


def test_score_variant_uses_sdk_and_caches(tmp_path, monkeypatch):
    pytest.importorskip("alphagenome")
    from alphagenome.models.v2 import variant_scorers

    tidy = pd.DataFrame(
        {
            "output_type": ["RNA_SEQ"],
            "track_name": ["t"],
            "gene_name": ["TP53"],
            "biosample_name": ["liver"],
            "raw_score": [0.5],
            "variant_scorer": [
                "GeneMaskLFCScorer(requested_output=RNA_SEQ, aggregation_type=DIFF_MEAN)"
            ],
        }
    )
    monkeypatch.setattr(variant_scorers, "tidy_scores", lambda scores: tidy)

    class FakeModel:
        calls = []

        def score_variant(self, interval, variant, condition):
            self.calls.append((interval.start, interval.end, variant.position, condition))
            return ["anndata placeholder"]

    client = AlphaGenomeClient("key", tmp_path, sequence_length=16_384, min_interval_s=0)
    client._model = FakeModel()
    variant = GenomicVariant.parse("chr17:7675088:C>T")

    first = client.score_variant(variant, ["UBERON:0002107"])
    second = client.score_variant(variant, "UBERON:0002107")
    assert FakeModel.calls == [(7675088 - 8192, 7675088 + 8192, 7675088, "UBERON:0002107")]
    assert not first.from_cache and second.from_cache
    assert second.summary[0].scorer == "GeneMaskLFCScorer(DIFF_MEAN)"
    assert second.label == "liver"


def _track(values, names=("t0", "t1")):
    from types import SimpleNamespace

    import numpy as np

    return SimpleNamespace(
        values=np.asarray(values, float), resolution=1, metadata=pd.DataFrame({"name": list(names)})
    )


def test_raw_track_predictions(tmp_path):
    pytest.importorskip("alphagenome")
    from types import SimpleNamespace

    from geneme_structure.models import GenomicInterval

    class FakeModel:
        def predict_variant(self, interval, variant, condition, requested_outputs):
            ref = SimpleNamespace(rna_seq=_track([[1, 1], [1, 1]]), dnase=None)
            alt = SimpleNamespace(rna_seq=_track([[1, 1], [1, 5]]), dnase=None)
            return SimpleNamespace(reference=ref, alternate=alt)

        def predict_interval(self, interval, condition, requested_outputs):
            return SimpleNamespace(rna_seq=_track([[0, 2], [0, 4]]), dnase=None)

    client = AlphaGenomeClient("key", tmp_path, sequence_length=16_384, min_interval_s=0)
    client._model = FakeModel()

    pred = client.predict_variant(
        GenomicVariant.parse("chr1:50000:A>G"), "UBERON:1", ["RNA_SEQ", "DNASE"]
    )
    assert [t.output_type for t in pred.tracks] == ["RNA_SEQ"]  # DNASE missing -> skipped
    assert (pred.tracks[0].top_track, pred.tracks[0].max_abs_delta) == ("t1", 4.0)

    region = GenomicInterval.parse("chr1:40001-60000")
    interval_pred = client.predict_interval(region, "UBERON:1", ["RNA_SEQ"])
    assert interval_pred.interval.width == 16_384
    assert (interval_pred.tracks[0].mean_signal, interval_pred.tracks[0].top_track) == (1.5, "t1")

    with pytest.raises(ValueError):
        client.predict_variant(GenomicVariant.parse("chr1:50000:A>G"), "UBERON:1", ["NOPE"])


def test_run_region_and_summary(settings, capsys):
    from geneme_structure.models import GenomicInterval
    from geneme_structure.pipeline import print_summary

    session = _session()
    session.add(
        "/overlap/region/",
        FakeResponse(
            payload=[
                {"id": "E1", "external_name": "TP53", "biotype": "protein_coding"},
                {"id": "E2", "external_name": "ORPHAN", "biotype": "protein_coding"},
            ]
        ),
    )
    session.add(
        "/uniprotkb/search",
        lambda m, u, kw: FakeResponse(
            payload={"results": [{"primaryAccession": "P04637"}]}
            if "TP53" in kw["params"]["query"]
            else {"results": []}
        ),
    )
    pipeline = GenomeToStructurePipeline(settings, use_alphagenome=False)
    pipeline.mapper = gene_mapping.GeneMapper(session=session)
    pipeline.alphafold = alphafold_client.AlphaFoldClient(settings.data_dir, session=session)

    result = pipeline.run_region(GenomicInterval.parse("chr17:7668421-7687490"), max_genes=5)
    assert result.genes_in_region == ["TP53", "ORPHAN"]
    assert [s.entry.entry_id for s in result.structures] == ["AF-P04637-F1"]
    assert result.warnings == ["No reviewed UniProt entry for ORPHAN"]

    print_summary(result)
    out = capsys.readouterr().out
    assert "Genes:       TP53, ORPHAN" in out and "! No reviewed UniProt entry" in out


def test_cli_gene_without_alphagenome(tmp_path, monkeypatch, capsys):
    session = _session()
    session.add("eutils", eutils([clinvar_doc("7", "NC_000017.11:7675087:C:T", protein="R2H")]))
    monkeypatch.setattr(gene_mapping, "build_session", lambda: session)
    monkeypatch.setattr(alphafold_client, "build_session", lambda: session)
    monkeypatch.setattr(clinvar.time, "sleep", lambda s: None)
    monkeypatch.chdir(tmp_path)

    code = main(["-o", str(tmp_path / "res"), "gene", "TP53", "--no-alphagenome"])
    assert code == 0
    out = capsys.readouterr().out
    assert "R2H" in out and "chr17:7675088:C>T" in out
    vcf = (tmp_path / "res" / "TP53_clinvar_missense.vcf").read_text()
    assert "chr17\t7675088\tR2H\tC\tT" in vcf
    summary = pd.read_csv(tmp_path / "res" / "TP53_clinvar_missense_summary.csv")
    assert summary.loc[0, "vcf_id"] == "R2H"
    prov = json.loads((tmp_path / "res" / "TP53_clinvar_missense_provenance.json").read_text())
    assert prov["data_sources"]["ensembl_release"] == 116
    assert prov["data_sources"]["clinvar_query_date"]
    viewer = (tmp_path / "res" / "TP53_clinvar_missense_viewer.html").read_text()
    assert "R2H" in viewer and "About this snapshot" in viewer


def test_cli_gene_no_hits_is_an_error(tmp_path, monkeypatch):
    session = _session().add("eutils", eutils([]))
    monkeypatch.setattr(gene_mapping, "build_session", lambda: session)
    monkeypatch.chdir(tmp_path)
    assert main(["-o", str(tmp_path), "gene", "NOTAGENE", "--no-alphagenome"]) == 2


def test_annotate_adds_allele_specific_clinvar(settings):
    vep = {
        **VEP_TP53,
        "colocated_variants": [
            {
                "id": "rs1",
                "clin_sig_allele": "T:pathogenic;T:benign",
                "var_synonyms": {"ClinVar": ["VCV000000007", "RCV1"]},
            }
        ],
    }
    session = FakeSession().add("/vep/", vep_echo(lambda line: vep))
    session.add(
        "eutils", eutils([clinvar_doc("7", "NC_000017.11:7675087:C:T", classification="Benign")])
    )
    pipeline = GenomeToStructurePipeline(settings, use_alphagenome=False)
    pipeline.mapper = gene_mapping.GeneMapper(session=session)
    pipeline.clinvar = clinvar.ClinVarClient(session=session)
    pipeline.clinvar.min_interval_s = 0
    [a] = pipeline.annotate([GenomicVariant.parse("chr17:7675088:C>T")])
    assert a.clinvar_label == "Benign" and a.clinvar_accession == "VCV000000007"

    broken = FakeSession().add("/vep/", vep_echo(lambda line: vep))
    broken.add("eutils", FakeResponse(status_code=500))
    pipeline.mapper = gene_mapping.GeneMapper(session=broken)
    pipeline.clinvar = clinvar.ClinVarClient(session=broken)
    pipeline.clinvar.min_interval_s = 0
    [a] = pipeline.annotate([GenomicVariant.parse("chr17:7675088:C>T")])
    assert a.clinvar_label == "pathogenic, benign (submissions)"  # falls back to VEP terms
