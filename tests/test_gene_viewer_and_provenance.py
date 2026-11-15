"""Gene-level 3D page and provenance metadata (offline)."""

import json

import pytest
from fakes import AF_ENTRY, TOY_PDB, FakeResponse, FakeSession

from geneme_structure import provenance as prov_mod
from geneme_structure.alphagenome_client import summarize_scores
from geneme_structure.gene_viewer import render_gene_html, variant_rows, write_gene_html
from geneme_structure.models import (
    AlphaFoldEntry,
    GenomicVariant,
    PipelineResult,
    ScoreRecord,
    StructureResult,
    VariantAnnotation,
    VariantScores,
)
from geneme_structure.provenance import collect_provenance, ensembl_release


def _results(tmp_path):
    pdb = tmp_path / "AF-P04637-F1-model_v6.pdb"
    pdb.write_text(TOY_PDB)
    entry = AlphaFoldEntry.model_validate({**AF_ENTRY, "latestVersion": 6})
    other = AlphaFoldEntry.model_validate({**AF_ENTRY, "entryId": "AF-P04637-F2"})
    variant = GenomicVariant.parse("chr17:7675088:C>T")
    recs = [
        ScoreRecord(
            output_type="RNA_SEQ",
            variant_scorer="GeneMaskLFCScorer(DIFF_MEAN)",
            gene_name="TP53",
            biosample_name="colon",
            raw_score=-0.25,
        ),
        ScoreRecord(
            output_type="SPLICE_SITES",
            variant_scorer="GeneMaskSplicingScorer(X)",
            gene_name="TP53",
            biosample_name="colon",
            raw_score=1.5,
        ),
    ]
    vs = VariantScores(
        query=str(variant),
        interval=variant.reference_interval.resized(16_384),
        condition="UBERON:0001157",
        summary=summarize_scores(recs),
        records=recs,
    )
    ann = VariantAnnotation(
        gene_symbol="TP53",
        consequence_terms=["missense_variant", "splice_region_variant"],
        clinvar_significance=["pathogenic"],
    )
    on = StructureResult(
        entry=entry,
        structure_path=str(pdb),
        mean_plddt=69.17,
        residue=2,
        residue_plddt=72.5,
        residue_confidence="confident",
    )
    off = on.model_copy(update={"entry": other})
    return [
        PipelineResult(
            query=str(variant), vcf_id="<R2H>", annotation=ann, structures=[on], variant_scores=[vs]
        ),
        PipelineResult(query="chr17:1:A>G", vcf_id="X9Y", structures=[off]),
        PipelineResult(query="chr17:2:A>G"),
    ]


def test_variant_rows(tmp_path):
    results = _results(tmp_path)
    rows = variant_rows(results, results[0].structures[0])
    assert rows[0]["consequence"] == "missense, splice region"
    assert rows[0]["residue"] == 2 and rows[0]["splicing"] == 1.5 and rows[0]["rna_lfc"] == -0.25
    assert rows[1]["residue"] is None  # other AlphaFold fragment: not on this model
    assert rows[2]["rna_lfc"] is None and rows[2]["label"] == "chr17:2:A>G"


def test_render_gene_html_escapes_and_embeds(tmp_path):
    results = _results(tmp_path)
    prov = {
        "generated_utc": "2026-01-01T00:00:00Z",
        "attribution": {"X": "<b>y</b>"},
        "data_sources": {
            "ensembl_release": 116,
            "alphafold_entries": {"AF-P04637-F1": {"model_version": 6}},
        },
    }
    page = render_gene_html(results, title="TP53 <test>", provenance=prov)
    assert "<title>TP53 &lt;test&gt;</title>" in page
    assert "&lt;R2H&gt;" in page and "<R2H>" not in page.split("<script>")[0]
    assert "AF-P04637-F1 v6" in page and "Ensembl release</dt><dd>116" in page
    assert "&lt;b&gt;y&lt;/b&gt;" in page
    assert "Mean pLDDT</dt><dd>69.2" in page
    assert "</script>" not in page.split("<script>")[-1].split("</script>")[0]
    out = write_gene_html(results, tmp_path / "x" / "v.html", title="t", provenance=prov)
    assert out.exists()


def test_render_gene_html_needs_a_structure(tmp_path):
    with pytest.raises(ValueError):
        render_gene_html([PipelineResult(query="q")], title="t", provenance={})


def test_collect_provenance(tmp_path, monkeypatch):
    session = FakeSession().add("/info/software", FakeResponse(payload={"release": 116}))
    monkeypatch.setattr(prov_mod, "git_commit", lambda: "abc1234")
    p = collect_provenance(
        _results(tmp_path), session=session, clinvar_date="2026-11-15", parameters={"gene": "TP53"}
    )
    json.dumps(p)  # serialisable
    ds = p["data_sources"]
    assert ds["ensembl_release"] == 116 and ds["clinvar_query_date"] == "2026-11-15"
    assert ds["alphafold_entries"]["AF-P04637-F1"]["model_version"] == 6
    assert ds["alphagenome"]["sequence_lengths"] == [16_384]
    assert p["software"]["git_commit"] == "abc1234" and p["parameters"] == {"gene": "TP53"}
    assert "AlphaGenome" in p["attribution"]


def test_ensembl_release_is_best_effort():
    down = FakeSession().add("/info/software", FakeResponse(status_code=503))
    assert ensembl_release(down) is None
    odd = FakeSession().add("/info/software", FakeResponse(payload={}))
    assert ensembl_release(odd) is None
