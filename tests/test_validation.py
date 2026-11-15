import math

import numpy as np
import pandas as pd
import pytest

from var2fold.clinvar import _tier, build_query
from var2fold.validation import (
    Feature,
    auc_by_group,
    bootstrap_auc_ci,
    clinvar_label,
    evaluate_features,
    mann_whitney_p,
    roc_auc,
    roc_points,
    vcf_clinvar_labels,
)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("Pathogenic", 1),
        ("Likely pathogenic", 1),
        ("Pathogenic/Likely pathogenic", 1),
        ("Likely_pathogenic", 1),
        ("Benign", 0),
        ("Benign/Likely benign", 0),
        ("Uncertain significance", None),
        ("Conflicting classifications of pathogenicity", None),
        ("", None),
        (None, None),
    ],
)
def test_clinvar_label(text, expected):
    assert clinvar_label(text) == expected


def test_roc_auc_perfect_chance_and_ties():
    assert roc_auc([0.9, 0.8, 0.2, 0.1], [1, 1, 0, 0]) == 1.0
    assert roc_auc([0.1, 0.2, 0.8, 0.9], [1, 1, 0, 0]) == 0.0
    assert roc_auc([0.5, 0.5, 0.5, 0.5], [1, 0, 1, 0]) == 0.5
    # One positive tied with one negative: 3 wins + 1 half out of 4 pairs.
    assert roc_auc([0.9, 0.5, 0.5, 0.1], [1, 1, 0, 0]) == pytest.approx(0.875)


def test_roc_auc_drops_nan_and_needs_both_classes():
    assert roc_auc([0.9, float("nan"), 0.1], [1, 1, 0]) == 1.0
    assert math.isnan(roc_auc([0.9, 0.8], [1, 1]))
    with pytest.raises(ValueError):
        roc_auc([0.1, 0.2], [1])


def test_roc_auc_matches_pair_counting():
    rng = np.random.default_rng(1)
    s = rng.integers(0, 6, 40).astype(float)  # many ties
    y = rng.integers(0, 2, 40)
    pos, neg = s[y == 1], s[y == 0]
    pairs = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    assert roc_auc(s, y) == pytest.approx(pairs / (len(pos) * len(neg)))


def test_mann_whitney_p():
    strong = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1] * 3
    labels = [1] * 5 + [0] * 5
    assert mann_whitney_p(strong, labels * 3) < 0.001
    assert mann_whitney_p([1, 1, 1, 1], [1, 0, 1, 0]) == 1.0
    assert mann_whitney_p([1, 2, 3, 4], [1, 0, 0, 1]) > 0.5


def test_bootstrap_ci_brackets_auc_and_is_deterministic():
    rng = np.random.default_rng(2)
    y = np.r_[np.ones(30, int), np.zeros(30, int)]
    s = np.r_[rng.normal(1, 1, 30), rng.normal(0, 1, 30)]
    lo, hi = bootstrap_auc_ci(s, y, n_boot=500, seed=3)
    assert lo < roc_auc(s, y) < hi
    assert (lo, hi) == bootstrap_auc_ci(s, y, n_boot=500, seed=3)
    assert all(math.isnan(x) for x in bootstrap_auc_ci([1.0], [1]))


def test_roc_points_endpoints():
    fpr, tpr = roc_points([0.9, 0.8, 0.2, 0.1], [1, 0, 1, 0])
    assert (fpr[0], tpr[0]) == (0.0, 0.0)
    assert (fpr[-1], tpr[-1]) == (1.0, 1.0)
    assert np.all(np.diff(fpr) >= 0) and np.all(np.diff(tpr) >= 0)


def _table():
    return pd.DataFrame(
        {
            "gene": ["A"] * 4 + ["B"] * 4,
            "label": [1, 1, 0, 0] * 2,
            "residue_plddt": [95, 90, 40, 30, 92, 50, 60, 55],
            "splicing": [-0.5, 0.4, 0.01, -0.02, 0.0, 0.0, 0.0, 0.0],
            "gnomad_af": [None, None, 0.01, 0.2, None, 1e-5, 0.3, None],
        }
    )


def test_feature_orientation_absolute_and_fill():
    t = _table()
    af = Feature("gnomad_af", "AF", higher_is_pathogenic=False, fill_missing=0.0)
    spl = Feature("splicing", "splicing", absolute=True)
    out = evaluate_features(t, [Feature("residue_plddt", "pLDDT"), spl, af], n_boot=200)
    by = out.set_index("feature")
    assert by.loc["residue_plddt", "auc"] == pytest.approx(0.875)
    # |-0.5| counts as large; gene B's all-zero scores tie (10 of 16 pairs won)
    assert by.loc["splicing", "auc"] == pytest.approx(0.625)
    # missing AF -> 0 -> rare, the pathogenic side (13.5 of 16 pairs)
    assert by.loc["gnomad_af", "auc"] == pytest.approx(0.84375)
    assert set(by["n_pathogenic"]) == {4} and set(by["n_benign"]) == {4}


def test_auc_by_group():
    t = _table()
    g = auc_by_group(t, features=[Feature("residue_plddt", "pLDDT")])
    assert g.loc["A", "residue_plddt"] == 1.0
    assert g.loc["B", "residue_plddt"] == 0.5


def test_vcf_clinvar_labels(tmp_path):
    vcf = tmp_path / "b.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7675088\tR175H\tC\tT\t.\t.\t"
        "CLNVCV=VCV1;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n"
        "chr17\t7676154\tP72R\tG\tC\t.\t.\tCLNVCV=VCV2;CLNSIG=Benign;CLNREVSTAT=x\n"
    )
    t = vcf_clinvar_labels(vcf)
    assert list(t["label"]) == [1, 0]
    assert t.loc[0, "review_status"] == "reviewed by expert panel"
    assert t.loc[0, "variant"].startswith("chr17:7675088")


def test_benign_significance_and_tier():
    q = build_query("TP53", "missense", ["benign", "likely_benign"])
    assert "clinsig_benign[prop] OR clinsig_likely_benign[prop]" in q
    firm = {"review_status": "reviewed by expert panel", "classification": "Benign"}
    likely = {"review_status": "reviewed by expert panel", "classification": "Likely benign"}
    assert _tier(firm) < _tier(likely)
