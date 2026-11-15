"""Plot layout helpers, UniProt feature parsing and data-driven findings."""

import pandas as pd
import pytest

from var2fold.gene_mapping import parse_uniprot_features
from var2fold.models import ProteinFeature
from var2fold.pipeline import batch_findings
from var2fold.plots import assign_rows, domain_lanes, group_labels


def test_group_labels_merges_shared_residue():
    assert group_labels({"R248Q": 248, "R248W": 248, "R175H": 175}) == {
        175: "R175H",
        248: "R248Q/W",
    }
    assert group_labels({"A1B": 5, "VCV1": 5}) == {5: "A1B/VCV1"}


def test_assign_rows_avoids_overlap_and_clamps():
    items = [(10.0, "AAAA"), (11.0, "BBBB"), (12.0, "CCCC"), (50.0, "DDDD")]
    placed = assign_rows(items, char_width=1.0, n_rows=2, bounds=(0.0, 51.0))
    assert [p[0] if p else None for p in placed] == [0, 1, None, 0]
    assert placed[3][1] == 49.0  # right edge 52 would cross 51 -> shifted left
    clamped = assign_rows([(0.0, "XXXX")], 1.0, 1, bounds=(0.0, 100.0))
    assert clamped[0] == (0, 2.0)


def test_domain_lanes_pack_overlaps():
    feats = [
        ProteinFeature(start=1, end=50, name="a", kind="Domain"),
        ProteinFeature(start=40, end=80, name="b", kind="Domain"),
        ProteinFeature(start=60, end=90, name="c", kind="Domain"),
    ]
    assert domain_lanes(feats) == [0, 1, 0]


def test_parse_uniprot_features_keeps_architecture():
    entry = {
        "features": [
            {
                "type": "Region",
                "description": "Interaction with X",
                "location": {"start": {"value": 1}, "end": {"value": 9}},
            },
            {
                "type": "Region",
                "description": "Disordered",
                "location": {"start": {"value": 1}, "end": {"value": 9}},
            },
            {
                "type": "Region",
                "description": "Oligomerization",
                "location": {"start": {"value": 325}, "end": {"value": 356}},
            },
            {
                "type": "DNA binding",
                "description": "",
                "location": {"start": {"value": 102}, "end": {"value": 292}},
            },
            {
                "type": "Domain",
                "description": "Kinase",
                "location": {"start": {"value": 5}, "end": {"value": None}},
            },
        ]
    }
    feats = parse_uniprot_features(entry)
    assert [(f.name, f.start) for f in feats] == [("DNA binding", 102), ("Oligomerization", 325)]


def _table():
    return pd.DataFrame(
        {
            "variant": ["v1", "v2", "v3", "v4"],
            "vcf_id": ["A1B", "C2D", "E3F", None],
            "residue_plddt": [95.0, 80.0, 40.0, None],
            "condition": ["colon"] * 4,
            "rna_lfc": [0.01, -0.02, 0.005, 0.0],
            "splicing": [0.5, 0.1, 0.05, 0.02],
            "consequence": [
                "missense_variant&splice_region_variant",
                "missense_variant",
                None,
                None,
            ],
        }
    )


def test_batch_findings_describes_data():
    lines = batch_findings(_table())
    assert lines[0] == "4 variant(s) analysed."
    assert "1 very high (>90), 1 confident (70-90), 1 low or very low (<=70)" in lines[1]
    assert "Strongest expression effect (RNA-seq gene LFC): C2D (-0.020)" in lines[2]
    assert "more than twice the next (C2D, +0.100)" in lines[3]
    assert "A1B is annotated by VEP as splice-region" in lines[3]
    assert batch_findings(pd.DataFrame()) == ["No variants were analysed."]


def test_figures_render():
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    from var2fold.plots import effects_panels, plddt_profile

    plddt = {i: 90.0 if 20 < i < 60 else 40.0 for i in range(1, 101)}
    feats = [ProteinFeature(start=20, end=60, name="Core", kind="Domain")]
    fig = plddt_profile(plddt, {"A30B": 30, "C31D": 31, "E80F": 80}, feats, title="t")
    assert fig.axes[0].get_title(loc="left") == "t"
    fig2 = effects_panels(_table().assign(label=lambda t: t["variant"]))
    assert len(fig2.axes) == 2  # accessibility column absent -> no panel
    with pytest.raises(ValueError):
        effects_panels(pd.DataFrame({"label": ["x"]}))
