"""Benchmark helpers: do var2fold's features separate pathogenic from benign variants?

Pure functions over a labelled `batch_table` (no network, no extra dependencies):

* `clinvar_label` turns a ClinVar classification into 1 (pathogenic side),
  0 (benign side) or None (VUS, conflicting, anything else).
* `roc_auc` is the Mann-Whitney formulation (ties count one half), so it needs
  no scikit-learn; `roc_points` gives the curve for plotting.
* `bootstrap_auc_ci` resamples within each class, keeping both classes present.
* `evaluate_features` builds the summary table: one row per feature, with AUC,
  95% bootstrap CI and a two-sided Mann-Whitney p-value (normal approximation).

AUC is oriented by each feature's stated direction ("higher means pathogenic"
or "lower means pathogenic"), so 0.5 is chance and 1.0 is perfect separation.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Feature:
    column: str
    description: str
    higher_is_pathogenic: bool = True
    absolute: bool = False  # score the magnitude, ignoring the sign of the effect
    fill_missing: float | None = None  # e.g. absent from gnomAD -> allele frequency 0

    def values(self, table: pd.DataFrame) -> pd.Series:
        v = pd.to_numeric(table[self.column], errors="coerce")
        if self.fill_missing is not None:
            v = v.fillna(self.fill_missing)
        return v.abs() if self.absolute else v


DEFAULT_FEATURES: tuple[Feature, ...] = (
    Feature("residue_plddt", "AlphaFold pLDDT at the mutated residue"),
    Feature("splicing", "AlphaGenome |splicing effect|", absolute=True),
    Feature("rna_lfc", "AlphaGenome |expression effect| (RNA-seq gene LFC)", absolute=True),
    Feature("accessibility", "AlphaGenome |accessibility effect| (DNase/ATAC)", absolute=True),
    Feature(
        "gnomad_af",
        "gnomAD allele frequency (reference only: ClinVar uses it to classify)",
        higher_is_pathogenic=False,
        fill_missing=0.0,
    ),
)


def clinvar_label(classification: str | None) -> int | None:
    """1 for Pathogenic / Likely pathogenic, 0 for Benign / Likely benign, else None."""
    if not classification:
        return None
    text = classification.lower().replace("_", " ")
    if "conflicting" in text or "uncertain" in text:
        return None
    pathogenic, benign = "pathogenic" in text, "benign" in text
    if pathogenic == benign:  # neither, or both (should not happen for aggregate calls)
        return None
    return 1 if pathogenic else 0


def _clean(scores: Sequence[float], labels: Sequence[int]) -> tuple[np.ndarray, np.ndarray]:
    s = np.asarray(scores, dtype=float)
    y = np.asarray(labels, dtype=int)
    if s.shape != y.shape:
        raise ValueError("scores and labels must have the same length")
    keep = ~np.isnan(s)
    return s[keep], y[keep]


def _rank_sum_u(s: np.ndarray, y: np.ndarray) -> tuple[float, int, int]:
    """Mann-Whitney U for the positive class, with average ranks for ties."""
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    sorted_s = s[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    n_pos, n_neg = int(y.sum()), int(len(y) - y.sum())
    u = ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2
    return float(u), n_pos, n_neg


def roc_auc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """P(score of a random positive > score of a random negative); ties count 1/2.

    NaN scores are dropped. Returns NaN unless both classes are present.
    """
    s, y = _clean(scores, labels)
    u, n_pos, n_neg = _rank_sum_u(s, y)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return u / (n_pos * n_neg)


def mann_whitney_p(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Two-sided p-value for U, normal approximation with tie correction."""
    s, y = _clean(scores, labels)
    u, n1, n2 = _rank_sum_u(s, y)
    if n1 == 0 or n2 == 0:
        return float("nan")
    n = n1 + n2
    _, counts = np.unique(s, return_counts=True)
    tie_term = float(((counts**3) - counts).sum()) / (n * (n - 1)) if n > 1 else 0.0
    var = n1 * n2 / 12 * ((n + 1) - tie_term)
    if var <= 0:
        return 1.0
    z = (abs(u - n1 * n2 / 2) - 0.5) / math.sqrt(var)  # continuity correction
    return float(min(1.0, math.erfc(max(z, 0.0) / math.sqrt(2))))


def bootstrap_auc_ci(
    scores: Sequence[float],
    labels: Sequence[int],
    n_boot: int = 2000,
    level: float = 0.95,
    seed: int = 0,
) -> tuple[float, float]:
    """Percentile CI for the AUC, resampling positives and negatives separately."""
    s, y = _clean(scores, labels)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    aucs = np.empty(n_boot)
    lab = np.r_[np.ones(len(pos), int), np.zeros(len(neg), int)]
    for b in range(n_boot):
        sample = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        aucs[b] = roc_auc(sample, lab)
    alpha = (1 - level) / 2
    return float(np.quantile(aucs, alpha)), float(np.quantile(aucs, 1 - alpha))


def roc_points(scores: Sequence[float], labels: Sequence[int]) -> tuple[np.ndarray, np.ndarray]:
    """(false positive rate, true positive rate) for every threshold, from (0,0) to (1,1)."""
    s, y = _clean(scores, labels)
    n_pos, n_neg = int(y.sum()), int(len(y) - y.sum())
    if n_pos == 0 or n_neg == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])
    thresholds = np.unique(s)[::-1]
    fpr = [0.0] + [float(((s >= t) & (y == 0)).sum()) / n_neg for t in thresholds]
    tpr = [0.0] + [float(((s >= t) & (y == 1)).sum()) / n_pos for t in thresholds]
    return np.array(fpr), np.array(tpr)


def oriented(feature: Feature, table: pd.DataFrame) -> pd.Series:
    """Feature values flipped where needed so that higher always means 'more pathogenic'."""
    v = feature.values(table)
    return v if feature.higher_is_pathogenic else -v


def evaluate_features(
    table: pd.DataFrame,
    features: Sequence[Feature] = DEFAULT_FEATURES,
    label_col: str = "label",
    n_boot: int = 2000,
    seed: int = 0,
) -> pd.DataFrame:
    """One row per feature: class sizes, AUC with 95% CI, and Mann-Whitney p."""
    rows = []
    y = table[label_col].astype(int)
    for f in features:
        s = oriented(f, table)
        keep = s.notna()
        lo, hi = bootstrap_auc_ci(s[keep], y[keep], n_boot=n_boot, seed=seed)
        rows.append(
            {
                "feature": f.column,
                "description": f.description,
                "n_pathogenic": int((y[keep] == 1).sum()),
                "n_benign": int((y[keep] == 0).sum()),
                "auc": roc_auc(s[keep], y[keep]),
                "ci_low": lo,
                "ci_high": hi,
                "p_value": mann_whitney_p(s[keep], y[keep]),
            }
        )
    return pd.DataFrame(rows)


def auc_by_group(
    table: pd.DataFrame,
    group_col: str = "gene",
    features: Sequence[Feature] = DEFAULT_FEATURES,
    label_col: str = "label",
) -> pd.DataFrame:
    """AUC per group (rows) and feature (columns); NaN where a class is missing."""
    out = {}
    for group, sub in table.groupby(group_col, sort=True):
        y = sub[label_col].astype(int)
        out[group] = {f.column: roc_auc(oriented(f, sub), y) for f in features}
    return pd.DataFrame.from_dict(out, orient="index")


def vcf_clinvar_labels(path) -> pd.DataFrame:
    """`variant`, `vcf_id`, `clnsig`, `review_status` and `label` for each row of a
    VCF written by `clinvar.write_vcf`; `variant` matches `batch_table`'s column."""
    from var2fold.models import GenomicVariant

    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom, pos, vid, ref, alt, _q, _f, info = line.rstrip("\n").split("\t")[:8]
            fields = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            clnsig = fields.get("CLNSIG", "").replace("_", " ")
            variant = GenomicVariant(
                chromosome=chrom, position=int(pos), reference_bases=ref, alternate_bases=alt
            )
            rows.append(
                {
                    "variant": str(variant),
                    "vcf_id": None if vid == "." else vid,
                    "clnsig": clnsig,
                    "review_status": fields.get("CLNREVSTAT", "").replace("_", " "),
                    "label": clinvar_label(clnsig),
                }
            )
    return pd.DataFrame(rows)
