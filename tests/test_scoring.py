import pandas as pd
import pytest

from var2fold.alphagenome_client import (
    records_from_tidy,
    scorer_kind,
    short_scorer_name,
    summarize_scores,
    variant_window,
)
from var2fold.models import GenomicInterval, GenomicVariant


def _tidy() -> pd.DataFrame:
    """Minimal frame with the columns `variant_scorers.tidy_scores` returns."""
    lfc = "GeneMaskLFCScorer(requested_output=RNA_SEQ, aggregation_type=DIFF_MEAN)"
    act = "GeneMaskActiveScorer(requested_output=RNA_SEQ)"
    atac = "CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=DIFF_LOG2_SUM)"
    return pd.DataFrame(
        {
            "variant_id": ["chr1:100:A>G"] * 5,
            "output_type": ["RNA_SEQ", "RNA_SEQ", "RNA_SEQ", "ATAC", "ATAC"],
            "variant_scorer": [lfc, lfc, act, atac, atac],
            "track_name": ["t1", "t2", "t1", "a1", None],
            "gene_name": ["TP53", "TP53", "TP53", None, None],
            "biosample_name": ["colon"] * 5,
            "raw_score": [0.1, -0.5, 3.0, 0.02, float("nan")],
            "extra_column": range(5),
        }
    )


def test_short_scorer_name_and_kind():
    long = "CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=ACTIVE_SUM)"
    assert short_scorer_name(long) == "CenterMaskScorer(width=501, ACTIVE_SUM)"
    assert short_scorer_name("GeneMaskActiveScorer(requested_output=RNA_SEQ)") == (
        "GeneMaskActiveScorer()"
    )
    assert scorer_kind(long) == "activity"
    assert scorer_kind("GeneMaskLFCScorer(DIFF_MEAN)") == "effect"


def test_records_drop_nan_scores_and_extra_columns():
    records = records_from_tidy(_tidy())
    assert len(records) == 4
    assert records[0].variant_scorer == "GeneMaskLFCScorer(DIFF_MEAN)"
    assert records[3].gene_name is None
    assert records_from_tidy(None) == [] and records_from_tidy(pd.DataFrame()) == []


def test_summarize_picks_max_abs_and_lists_effects_first():
    summary = summarize_scores(records_from_tidy(_tidy()))
    assert [s.kind for s in summary] == ["effect", "effect", "activity"]
    rna = next(s for s in summary if s.scorer.startswith("GeneMaskLFC"))
    assert (rna.n_scores, rna.top_score, rna.top_track, rna.top_gene) == (2, -0.5, "t2", "TP53")
    atac = next(s for s in summary if s.output_type == "ATAC")
    assert atac.top_gene is None


@pytest.mark.parametrize(
    ("variant", "expected_start"),
    [("chr17:7675088:C>T", 7675088 - 8192), ("chr1:1000000:ACG>A", 1000001 - 8192)],
)
def test_variant_window_matches_sdk_convention(variant, expected_start):
    # Values checked against alphagenome.data.genome.Variant.reference_interval.resize
    w = variant_window(GenomicVariant.parse(variant), 16_384)
    assert (w.start, w.width) == (expected_start, 16_384)


def test_resized_rejects_negative_start():
    with pytest.raises(ValueError):
        GenomicInterval(chromosome="chr1", start=10, end=11).resized(16_384)


def _scores(condition: str, biosample: str, rna_score: float, with_atac: bool = True):
    from var2fold.models import ScoreRecord, VariantScores

    lfc = "GeneMaskLFCScorer(DIFF_MEAN)"
    records = [
        ScoreRecord(
            output_type="RNA_SEQ",
            variant_scorer=lfc,
            track_name="t",
            gene_name="TP53",
            biosample_name=biosample,
            raw_score=rna_score,
        )
    ]
    if with_atac:
        records.append(
            ScoreRecord(
                output_type="ATAC",
                variant_scorer="CenterMaskScorer()",
                biosample_name=biosample,
                raw_score=0.2,
            )
        )
    return VariantScores(
        query="chr1:100:A>G",
        interval=GenomicInterval(chromosome="chr1", start=0, end=16_384),
        condition=condition,
        summary=summarize_scores(records),
        records=records,
    )


def test_compare_conditions_pivots_and_keeps_order():
    from var2fold.alphagenome_client import compare_conditions

    liver = _scores("UBERON:0002107", "liver", -0.3, with_atac=False)
    colon = _scores("UBERON:0001157", "transverse colon", 0.1)
    table = compare_conditions([liver, colon])
    assert list(table.columns) == ["liver", "transverse colon"]
    assert table.loc[("RNA_SEQ", "GeneMaskLFCScorer(DIFF_MEAN)"), "liver"] == -0.3
    assert pd.isna(table.loc[("ATAC", "CenterMaskScorer()"), "liver"])
    assert compare_conditions([]).empty


def test_label_falls_back_to_condition():
    vs = _scores("UBERON:0002107", "liver", 0.1)
    assert vs.label == "liver"
    assert vs.model_copy(update={"records": []}).label == "UBERON:0002107"
