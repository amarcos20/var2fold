import pytest
from pydantic import ValidationError

from geneme_structure.models import AlphaFoldEntry, GenomicInterval, GenomicVariant, plddt_band


def test_parse_variant_normalizes():
    v = GenomicVariant.parse("17-7675088-c-t")
    assert v.chromosome == "chr17"
    assert (v.position, v.reference_bases, v.alternate_bases) == (7675088, "C", "T")
    assert str(v) == "chr17:7675088:C>T"
    assert v.ensembl_chromosome == "17"


@pytest.mark.parametrize("bad", ["chr17:0:C>T", "chr99:10:C>T", "chr1:10:C>C", "chr1:10:C>Z"])
def test_invalid_variants(bad):
    with pytest.raises((ValueError, ValidationError)):
        GenomicVariant.parse(bad)


def test_interval_parse_is_one_based_inclusive():
    iv = GenomicInterval.parse("chrX:1,001-2,000")
    assert (iv.start, iv.end, iv.width) == (1000, 2000, 1000)
    assert str(iv) == "chrX:1001-2000"


def test_alphafold_entry_ignores_unknown_fields():
    e = AlphaFoldEntry.model_validate(
        {
            "entryId": "AF-P04637-F1",
            "uniprotAccession": "P04637",
            "uniprotStart": 1,
            "uniprotEnd": 393,
            "brandNewField": 42,
        }
    )
    assert e.covers(175) and not e.covers(400)


def test_plddt_band():
    assert [plddt_band(x) for x in (95, 80, 60, 10)] == [
        "very high",
        "confident",
        "low",
        "very low",
    ]
