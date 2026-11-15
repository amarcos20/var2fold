"""Offline tests for the REST clients, using tests/fakes.py."""

from pathlib import Path

import pytest
import requests
from fakes import AF_ENTRY, TOY_PDB, FakeResponse, vep_echo

from var2fold.alphafold_client import AlphaFoldClient
from var2fold.gene_mapping import VEP_BATCH_SIZE, GeneMapper
from var2fold.models import GenomicInterval, GenomicVariant

# ------------------------------------------------------------------ AlphaFold DB


def _af_client(tmp_path: Path, session) -> AlphaFoldClient:
    return AlphaFoldClient(tmp_path, session=session)


def test_get_predictions_parses_and_memoises(tmp_path, fake_session):
    fake_session.add("/api/prediction/P04637", FakeResponse(payload=[AF_ENTRY]))
    client = _af_client(tmp_path, fake_session)
    entries = client.get_predictions("p04637")
    assert entries[0].entry_id == "AF-P04637-F1" and entries[0].global_plddt == 69.2
    client.get_predictions("P04637")
    assert len(fake_session.calls) == 1  # second call served from memo


def test_get_predictions_404_is_empty(tmp_path, fake_session):
    fake_session.add("/api/prediction/", FakeResponse(status_code=404))
    assert _af_client(tmp_path, fake_session).get_predictions("Q9H9K5") == []


def test_get_predictions_rejects_bad_accession(tmp_path, fake_session):
    with pytest.raises(ValueError):
        _af_client(tmp_path, fake_session).get_predictions("not-an-id")
    assert fake_session.calls == []


def test_server_error_propagates(tmp_path, fake_session):
    fake_session.add("/api/prediction/", FakeResponse(status_code=500))
    with pytest.raises(requests.HTTPError):
        _af_client(tmp_path, fake_session).get_predictions("P04637")


def test_fetch_structure_downloads_once_and_reads_residue(tmp_path, fake_session):
    fake_session.add("/api/prediction/", FakeResponse(payload=[AF_ENTRY]))
    fake_session.add("/files/", FakeResponse(content=TOY_PDB.encode()))
    client = _af_client(tmp_path, fake_session)

    result = client.fetch_structure("P04637", residue=2)
    assert Path(result.structure_path).read_text() == TOY_PDB
    assert (result.residue_plddt, result.residue_confidence) == (72.5, "confident")
    assert result.mean_plddt == pytest.approx((95 + 72.5 + 40) / 3, abs=0.01)

    client.fetch_structure("P04637", residue=3)
    downloads = [c for c in fake_session.calls if "/files/" in c[1]]
    assert len(downloads) == 1  # cached on disk


def test_fragment_offset_maps_uniprot_residue_to_model_numbering(tmp_path, fake_session):
    fragment = {
        **AF_ENTRY,
        "entryId": "AF-P04637-F2",
        "uniprotStart": 101,
        "uniprotEnd": 103,
        "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P04637-F2-model_v6.pdb",
    }
    fake_session.add("/api/prediction/", FakeResponse(payload=[AF_ENTRY, fragment]))
    fake_session.add("/files/", FakeResponse(content=TOY_PDB.encode()))
    result = _af_client(tmp_path, fake_session).fetch_structure("P04637", residue=102)
    assert result.entry.entry_id == "AF-P04637-F2"  # the fragment covering residue 102
    assert result.residue_plddt == 72.5  # model residue 2


# ------------------------------------------------------------------ Ensembl / UniProt


def test_annotate_variants_batches_and_keeps_order(fake_session):
    def record(line: str) -> dict:
        pos = int(line.split()[1])
        return {
            "allele_string": "A/G",
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{"gene_symbol": f"G{pos}", "canonical": 1}],
        }

    fake_session.add("/vep/homo_sapiens/region", vep_echo(record))
    variants = [GenomicVariant.parse(f"chr1:{i}:A>G") for i in range(1, VEP_BATCH_SIZE + 6)]
    annotations = GeneMapper(session=fake_session).annotate_variants(variants)

    assert len(fake_session.calls) == 2  # 200 + 5
    assert [len(c[2]["json"]["variants"]) for c in fake_session.calls] == [VEP_BATCH_SIZE, 5]
    assert fake_session.calls[0][2]["json"]["uniprot"] == 1
    assert [a.gene_symbol for a in annotations[:3]] == ["G1", "G2", "G3"]
    assert annotations[-1].gene_symbol == f"G{VEP_BATCH_SIZE + 5}"


def test_annotate_variant_missing_record_gives_empty_annotation(fake_session):
    fake_session.add("/vep/", FakeResponse(payload=[]))
    a = GeneMapper(session=fake_session).annotate_variant(GenomicVariant.parse("chr1:5:A>G"))
    assert a.gene_symbol is None and a.uniprot_accessions == []


def test_genes_in_region_filters_protein_coding(fake_session):
    fake_session.add(
        "/overlap/region/homo_sapiens/17:101-200",
        FakeResponse(
            payload=[
                {"id": "ENSG1", "external_name": "TP53", "biotype": "protein_coding"},
                {"id": "ENSG2", "external_name": "WRAP53", "biotype": "protein_coding"},
                {"id": "ENSG3", "external_name": "MIR123", "biotype": "miRNA"},
                {"id": "ENSG1", "external_name": "TP53", "biotype": "protein_coding"},
            ]
        ),
    )
    mapper = GeneMapper(session=fake_session)
    region = GenomicInterval(chromosome="chr17", start=100, end=200)
    assert mapper.genes_in_region(region) == ["TP53", "WRAP53"]
    assert fake_session.calls[0][2]["params"] == {"feature": "gene"}


def test_genes_in_region_rejects_huge_regions(fake_session):
    with pytest.raises(ValueError):
        GeneMapper(session=fake_session).genes_in_region(
            GenomicInterval(chromosome="chr1", start=0, end=10_000_000)
        )


def test_uniprot_for_gene(fake_session):
    fake_session.add(
        "/uniprotkb/search",
        lambda m, u, kw: FakeResponse(
            payload={"results": [{"primaryAccession": "P04637"}]}
            if "gene_exact:TP53" in kw["params"]["query"]
            else {"results": []}
        ),
    )
    mapper = GeneMapper(session=fake_session)
    assert mapper.uniprot_for_gene("TP53") == "P04637"
    assert mapper.uniprot_for_gene("NOPE") is None
    assert "reviewed:true" in fake_session.calls[0][2]["params"]["query"]
