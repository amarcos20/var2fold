"""ClinVar gene lookup: query building, SPDI conversion, ranking and export (offline)."""

import pytest
from fakes import FakeResponse, FakeSession, clinvar_doc, eutils

from var2fold.clinvar import (
    ClinVarClient,
    allele_key,
    build_query,
    matches_consequence,
    parse_summary,
    select_spread,
    spdi_key,
    spdi_to_vcf,
    write_vcf,
)
from var2fold.models import GenomicVariant, VariantAnnotation
from var2fold.vcf import read_vcf


def test_build_query_filters():
    q = build_query("BRCA1", "missense", ["pathogenic"])
    assert q.startswith("BRCA1[gene] AND (clinsig_pathogenic[prop]) AND single_gene[prop]")
    assert '"missense variant"[molecular consequence]' in q
    assert "molecular consequence" not in build_query("BRCA1", "any", ["pathogenic"])


@pytest.mark.parametrize(
    "args",
    [
        ("BRCA1", "nonsense-ish", ["pathogenic"]),
        ("BRCA1", "missense", ["benign"]),
        ("BRCA1", "missense", []),
        ("BRCA1) OR (x", "missense", ["pathogenic"]),
    ],
)
def test_build_query_rejects_bad_input(args):
    with pytest.raises(ValueError):
        build_query(*args)


@pytest.mark.parametrize(
    ("spdi", "expected"),
    [
        ("NC_000017.11:7674220:G:A", (7674221, "G", "A")),  # SNV, 0-based -> 1-based
        ("NC_000007.14:117509089:GG:G", (117509090, "GG", "G")),  # repeat deletion
        ("NC_000007.14:117592386:AGATTCTGAG:AG", (117592387, "AGATTCTGA", "A")),
        ("NC_000007.14:117509115:TT:TTT", (117509116, "T", "TT")),  # duplication
        ("NC_000007.14:117536630:GC:AA", (117536631, "GC", "AA")),  # MNV untouched
        ("NC_000007.14:117592506:TCAGAAC:", (117592507, "TCAGAAC", "")),  # needs anchor
    ],
)
def test_spdi_to_vcf(spdi, expected):
    assert spdi_to_vcf(spdi) == expected


def test_parse_summary_rejects_unusable_records():
    good = clinvar_doc("1", "NC_000017.11:7674220:G:A")
    assert parse_summary(good, "tp53")["chromosome"] == "chr17"
    assert parse_summary(good, "BRCA1") is None  # other gene
    compound = clinvar_doc("2", "NC_000017.11:1:G:A")
    compound["variation_set"] *= 2
    assert parse_summary(compound, "TP53") is None
    no_38 = clinvar_doc("3", "NC_000017.11:1:G:A")
    no_38["variation_set"][0]["variation_loc"].pop(0)
    assert parse_summary(no_38, "TP53") is None


def test_missense_filter_drops_stop_gains_and_frameshifts():
    assert matches_consequence({"protein_change": "R175H"}, "missense")
    assert not matches_consequence({"protein_change": "E1817*"}, "missense")
    assert not matches_consequence({"protein_change": "K684fs"}, "missense")
    assert matches_consequence({"protein_change": "E1817*"}, "coding")


def _p(pos, review="reviewed by expert panel", cls="Pathogenic"):
    return {"chromosome": "chr1", "position": pos, "review_status": review, "classification": cls}


def test_select_spread_prefers_trust_then_spreads_ties():
    single = "criteria provided, single submitter"
    parsed = [_p(pos) for pos in range(100, 110)] + [_p(1, single), _p(2, cls="Likely pathogenic")]
    picked = select_spread(parsed, 4)
    assert [p["position"] for p in picked] == [101, 103, 106, 108]  # evenly across the tier
    everything = select_spread(parsed, None)
    assert [p["position"] for p in everything][-2:] == [2, 1]  # likely-path, then 1 star
    assert select_spread(parsed, 0) == []


def test_variants_for_gene_resolves_indel_anchor():
    docs = [
        clinvar_doc("11", "NC_000017.11:7674220:G:A", protein="R248W"),
        clinvar_doc("12", "NC_000017.11:7675205:CT:", protein="V1L"),
        clinvar_doc("13", "NC_000017.11:7675300:C:T", protein="Q3*"),  # dropped: not missense
        clinvar_doc("14", "NC_000017.11:7675400:C:T", review="no assertion criteria provided"),
    ]
    session = FakeSession().add("eutils", eutils(docs))
    session.add("/sequence/region/", FakeResponse(text="A\n"))
    client = ClinVarClient(session=session, api_key="secret-key")
    client.min_interval_s = 0
    picked = client.variants_for_gene("TP53", max_variants=2)

    assert [str(v.variant) for v in picked] == ["chr17:7674221:G>A", "chr17:7675205:ACT>A"]
    assert picked[0].label == "R248W" and picked[0].stars == 3
    seq_call = [c for c in session.calls if "/sequence/region/" in c[1]][0]
    assert seq_call[1].endswith("/17:7675205..7675205:1")
    # key travels in the POST body, never in the URL or repr
    eutils_calls = [c for c in session.calls if "eutils" in c[1]]
    assert all("secret" not in url for _, url, _ in eutils_calls)
    assert eutils_calls[0][2]["data"]["api_key"] == "secret-key"
    assert "secret" not in repr(client)


def test_variants_for_gene_empty_search():
    session = FakeSession().add("eutils", eutils([]))
    client = ClinVarClient(session=session)
    assert client.variants_for_gene("TP53") == []


def test_write_vcf_round_trips(tmp_path):
    session = FakeSession().add("eutils", eutils([clinvar_doc("11", "NC_000017.11:7674220:G:A")]))
    client = ClinVarClient(session=session)
    client.min_interval_s = 0
    picked = client.variants_for_gene("TP53")
    path = write_vcf(picked, tmp_path / "sel.vcf")
    text = path.read_text()
    assert "CLNREVSTAT=reviewed_by_expert_panel" in text
    [record] = list(read_vcf(path))
    assert record.variant == picked[0].variant
    assert record.vcf_id == "R2H"


def test_allele_key_matches_vcf_and_spdi():
    assert allele_key("17", 7676154, "G", "C") == spdi_key("17", "NC_000017.11:7676153:G:C")
    # left-aligned VCF deletion vs canonical (repeat-expanded) SPDI
    assert allele_key("chr7", 117509089, "AG", "A") == spdi_key("7", "NC_000007.14:117509089:GG:G")
    assert allele_key("7", 117531086, "T", "TT") == spdi_key("7", "NC_000007.14:117531085:T:TT")
    assert allele_key("17", 1, "G", "C") != allele_key("17", 1, "G", "A")


def test_classify_matches_exact_allele_only():
    docs = [
        clinvar_doc("12351", "NC_000017.11:7676153:G:C", classification="Benign"),
        clinvar_doc(
            "864652",
            "NC_000017.11:7676153:G:A",
            classification="Conflicting classifications of pathogenicity",
        ),
    ]
    session = FakeSession().add("eutils", eutils(docs))
    client = ClinVarClient(session=session)
    client.min_interval_s = 0
    p72r = GenomicVariant.parse("chr17:7676154:G>C")
    other = GenomicVariant.parse("chr17:7676154:G>T")  # no ClinVar record of its own
    hits = client.classify([(p72r, ["VCV000012351", "VCV000864652"]), (other, ["VCV000864652"])])
    assert hits[0] == {
        "classification": "Benign",
        "review_status": "reviewed by expert panel",
        "accession": "VCV000012351",
    }
    assert hits[1] is None
    assert client.classify([(p72r, [])]) == [None]


def test_clinvar_label_prefers_aggregate():
    a = VariantAnnotation(clinvar_significance=["pathogenic", "benign"])
    assert a.clinvar_label == "pathogenic, benign (submissions)"
    assert a.model_copy(update={"clinvar_classification": "Benign"}).clinvar_label == "Benign"
    assert VariantAnnotation().clinvar_label is None
