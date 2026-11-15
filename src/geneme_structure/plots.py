"""Figures for notebooks and reports (needs the `notebook` extra: matplotlib).

* `plddt_profile`: AlphaFold confidence along the sequence, the analysed
  variants marked with leader-line labels, and the UniProt domain architecture.
* `effects_panels`: one lollipop panel per AlphaGenome scorer family.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

from geneme_structure.models import ProteinFeature

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
SERIES = "#2a78d6"
PLDDT_BANDS = [
    (90, 100, "#0053D6", "very high"),
    (70, 90, "#65CBF3", "confident"),
    (50, 70, "#FFDB13", "low"),
    (0, 50, "#FF7D45", "very low"),
]
EFFECT_PANELS = [
    ("rna_lfc", "Expression (RNA-seq gene LFC)"),
    ("splicing", "Splicing (strongest splice scorer)"),
    ("accessibility", "Accessibility (DNase/ATAC log2 ratio)"),
]
LABEL_ROWS = (103, 110, 117, 124)  # variant label rows, in a strip above the 0-100 scale
DOMAIN_Y, DOMAIN_LANE = 8.0, 9.0


def style(ax: Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def group_labels(variants: Mapping[str, int]) -> dict[int, str]:
    """Merge variants that share a residue: {'R248Q': 248, 'R248W': 248} -> {248: 'R248Q/W'}."""
    by_residue: dict[int, list[str]] = {}
    for label, residue in variants.items():
        by_residue.setdefault(residue, []).append(label)
    out = {}
    for residue, labels in sorted(by_residue.items()):
        prefix = labels[0][:-1]
        if len(labels) > 1 and all(lab[:-1] == prefix for lab in labels):
            out[residue] = prefix + "/".join(lab[-1] for lab in labels)
        else:
            out[residue] = "/".join(labels)
    return out


def assign_rows(
    items: Sequence[tuple[float, str]],
    char_width: float,
    n_rows: int,
    gap: float = 0.0,
    bounds: tuple[float, float] | None = None,
) -> list[tuple[int, float] | None]:
    """Place centred labels in the first row where they do not overlap (x in data units).

    Labels are shifted inside `bounds` when they would spill past an edge.
    Returns (row, label_x) per item, in input order, or None when no row has room.
    """
    row_end = [float("-inf")] * n_rows
    out: list[tuple[int, float] | None] = [None] * len(items)
    for i in sorted(range(len(items)), key=lambda i: items[i][0]):
        x, text = items[i]
        half = len(text) * char_width / 2
        if bounds is not None:
            x = min(max(x, bounds[0] + half), bounds[1] - half)
        for r in range(n_rows):
            if x - half > row_end[r] + gap:
                row_end[r] = x + half
                out[i] = (r, x)
                break
    return out


def domain_lanes(features: Sequence[ProteinFeature]) -> list[int]:
    """Greedy interval packing: overlapping features go to separate lanes."""
    lane_end: list[int] = []
    lanes = []
    for f in features:
        for i, end in enumerate(lane_end):
            if f.start > end:
                lane_end[i] = f.end
                lanes.append(i)
                break
        else:
            lane_end.append(f.end)
            lanes.append(len(lane_end) - 1)
    return lanes


def plddt_profile(
    plddt: Mapping[int, float],
    variants: Mapping[str, int],
    domains: Sequence[ProteinFeature] = (),
    title: str | None = None,
    ax: Axes | None = None,
) -> Figure:
    """pLDDT along the sequence with variants and domain architecture.

    Variant labels sit in a strip above the curve, with leader lines down to
    each point; labels that would overlap are dropped and counted in the x label.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4.4), dpi=150)
    fig = ax.figure
    length = max(plddt)
    bounds = (1.0, float(length))

    for lo, hi, colour, name in PLDDT_BANDS:
        ax.axhspan(lo, hi, color=colour, alpha=0.12, linewidth=0)
        ax.text(length * 1.01, (lo + hi) / 2, name, va="center", fontsize=8, color=MUTED)
    ax.plot(list(plddt), list(plddt.values()), color=SERIES, linewidth=2)

    # Width of one character in residues, from the axes size (8 pt and 7.5 pt text).
    axes_px = ax.get_position().width * fig.get_figwidth() * fig.dpi
    char_width = 0.62 * (8 * fig.dpi / 72) * length / axes_px

    points = [(r, n) for r, n in group_labels(variants).items() if r in plddt]
    placed = assign_rows(points, char_width, len(LABEL_ROWS), char_width, bounds)
    unlabelled = 0
    for (residue, name), slot in zip(points, placed, strict=True):
        y = plddt[residue]
        ax.plot(
            residue,
            y,
            "o",
            markersize=7,
            color=SERIES,
            markeredgecolor="white",
            markeredgewidth=2,
            zorder=3,
        )
        if slot is None:
            unlabelled += 1
            continue
        row, lx = slot
        ly = LABEL_ROWS[row]
        ax.plot([residue, residue], [y + 2, ly - 1], color=MUTED, linewidth=0.8, zorder=2)
        # Surface-coloured box hides leader lines of higher rows passing behind.
        ax.text(
            lx,
            ly,
            name,
            ha="center",
            va="bottom",
            fontsize=8,
            color=INK,
            zorder=4,
            bbox={"boxstyle": "square,pad=0.1", "facecolor": SURFACE, "edgecolor": "none"},
        )

    lanes = domain_lanes(domains)
    for lane in sorted(set(lanes)):
        in_lane = [f for f, ln in zip(domains, lanes, strict=True) if ln == lane]
        y = DOMAIN_Y + lane * DOMAIN_LANE
        labels = assign_rows(
            [((f.start + f.end) / 2, f.name) for f in in_lane],
            char_width * 0.95,
            1,
            char_width,
            bounds,
        )
        for f, slot in zip(in_lane, labels, strict=True):
            ax.plot([f.start, f.end], [y, y], color=MUTED, linewidth=2, solid_capstyle="butt")
            if slot is not None:
                ax.text(slot[1], y + 2, f.name, ha="center", va="bottom", fontsize=7.5, color=MUTED)

    ax.set_xlim(1, length)
    ax.set_ylim(0, LABEL_ROWS[-1] + 6)
    ax.set_yticks(range(0, 101, 20))
    ax.spines["left"].set_bounds(0, 100)
    xlabel = "Residue"
    if unlabelled:
        xlabel += f"   ·   {unlabelled} variant(s) not labelled for space (see the table)"
    ax.set_xlabel(xlabel, color=MUTED, fontsize=9)
    ax.set_ylabel("pLDDT", color=MUTED, fontsize=9)
    if title:
        ax.set_title(title, loc="left", fontsize=11, color=INK)
    style(ax)
    fig.tight_layout()
    return fig


def effects_panels(
    table: pd.DataFrame,
    label_col: str = "label",
    title: str | None = None,
    panels: Sequence[tuple[str, str]] = EFFECT_PANELS,
) -> Figure:
    """One lollipop panel per scorer family; compare within a panel, never across."""
    import matplotlib.pyplot as plt

    panels = [(c, t) for c, t in panels if c in table and table[c].notna().any()]
    if not panels:
        raise ValueError("table has no AlphaGenome score columns (run with AlphaGenome)")
    data = table.drop_duplicates(label_col).set_index(label_col)
    order = list(data.index)[::-1]  # first row at the top
    data = data.loc[order]

    fig, axes = plt.subplots(
        1,
        len(panels),
        figsize=(3.7 * len(panels), 0.9 + 0.36 * len(order)),
        dpi=150,
        sharey=True,
        squeeze=False,
    )
    for ax, (col, panel_title) in zip(axes[0], panels, strict=True):
        values = data[col].astype(float).fillna(0.0)
        ax.axvline(0, color=MUTED, linewidth=1)
        ax.hlines(range(len(values)), 0, values, color=SERIES, linewidth=2)
        ax.plot(
            values,
            range(len(values)),
            "o",
            markersize=7,
            color=SERIES,
            markeredgecolor="white",
            markeredgewidth=2,
        )
        top = values.abs().idxmax()
        side = 1 if values[top] >= 0 else -1  # label beyond the dot, away from its stem
        ax.annotate(
            f"{values[top]:+.3f}",
            (values[top], order.index(top)),
            xytext=(7 * side, 0),
            textcoords="offset points",
            va="center",
            ha="left" if side > 0 else "right",
            fontsize=8,
            color=INK,
        )
        ax.set_title(panel_title, loc="left", fontsize=10, color=INK)
        ax.set_yticks(range(len(order)), order)
        lim = (values.abs().max() or 1.0) * 1.5
        ax.set_xlim(-lim, lim)
        style(ax)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", color=GRID, linewidth=0.8)
    if title:
        fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    return fig
