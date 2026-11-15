"""One self-contained HTML page per gene: every analysed variant on the AlphaFold model.

Click a variant in the table to fly to its residue; the side panel shows its
annotation, confidence and AlphaGenome scores. A provenance panel records the
date, software and database versions, and the data-source terms.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from geneme_structure.models import PipelineResult, StructureResult
from geneme_structure.pipeline import _top_effect
from geneme_structure.report import PLDDT_BANDS, THREEDMOL_URL, _e, _script_json

_PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1d2330; --muted: #5d6676; --panel: #f4f6f9;
    --border: #dde2ea; --accent: #c2185b; --row: #fbe9f0; --viewer-bg: #ffffff;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #12151b; --fg: #e6e9ef; --muted: #9aa3b2; --panel: #1b2029;
      --border: #2c3340; --accent: #f06292; --row: #3a1d2a; --viewer-bg: #12151b;
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: var(--bg); color: var(--fg);
    font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  main {{ max-width: 1240px; margin: 0 auto; padding: 24px 16px 48px; }}
  nav {{ font-size: .85rem; margin-bottom: 12px; }}
  h1 {{ font-size: 1.5rem; margin: 0 0 4px; }}
  .sub {{ color: var(--muted); margin: 0 0 20px; }}
  .layout {{ display: grid; grid-template-columns: minmax(0, 2fr) minmax(0, 1fr); gap: 20px; }}
  @media (max-width: 860px) {{ .layout {{ grid-template-columns: 1fr; }} }}
  #viewer {{
    position: relative; width: 100%; height: 540px; border: 1px solid var(--border);
    border-radius: 10px; overflow: hidden; background: var(--viewer-bg);
  }}
  @media (max-width: 600px) {{ #viewer {{ height: 380px; }} }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }}
  button {{
    font: inherit; padding: 6px 12px; border-radius: 6px; cursor: pointer;
    border: 1px solid var(--border); background: var(--panel); color: var(--fg);
  }}
  button:hover, button:focus-visible {{ border-color: var(--accent); }}
  .panel {{
    background: var(--panel); border: 1px solid var(--border);
    border-radius: 10px; padding: 14px 16px; margin-bottom: 16px;
  }}
  .panel h2 {{
    font-size: .8rem; text-transform: uppercase; letter-spacing: .06em;
    color: var(--muted); margin: 0 0 8px;
  }}
  section .panel {{ margin-top: 16px; }}
  dl {{ display: grid; grid-template-columns: auto 1fr; gap: 4px 12px; margin: 0; }}
  dt {{ color: var(--muted); }}
  dd {{ margin: 0; overflow-wrap: anywhere; }}
  .legend span {{
    display: inline-block; width: 12px; height: 12px; border-radius: 2px;
    margin-right: 6px; vertical-align: -1px;
  }}
  .legend div {{ margin: 2px 0; }}
  .table-wrap {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; width: 100%; font-size: .85rem; }}
  th, td {{ text-align: left; padding: 5px 6px; border-bottom: 1px solid var(--border); }}
  th {{ color: var(--muted); font-weight: 600; }}
  td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  tbody tr {{ cursor: pointer; }}
  tbody tr:hover {{ background: var(--bg); }}
  tbody tr.sel {{ background: var(--row); }}
  td button {{ padding: 1px 8px; font-weight: 600; }}
  .off td {{ color: var(--muted); }}
  .note {{ color: var(--muted); font-size: .8rem; margin: 8px 0 0; }}
  .hl {{ color: var(--accent); font-weight: 600; }}
  .prov dl {{ font-size: .85rem; }}
  .prov ul {{ margin: 6px 0 0; padding-left: 18px; font-size: .8rem; color: var(--muted); }}
  .warn {{ border-left: 3px solid var(--accent); padding-left: 10px; font-size: .85rem; }}
  a {{ color: inherit; }}
</style>
</head>
<body>
<main>
  {nav}
  <h1>{heading}</h1>
  <p class="sub">{subheading}</p>
  <div class="layout">
    <section>
      <div id="viewer" role="img" aria-label="3D structure of {entry_id} with variants"></div>
      <div class="controls">
        <button id="btn-all">Whole protein</button>
        <button id="btn-labels" aria-pressed="false">Label all variants</button>
        <button id="btn-spin" aria-pressed="false">Spin</button>
        <button id="btn-surface" aria-pressed="false">Surface</button>
      </div>
      <div class="panel">
        <h2>Variants</h2>
        <div class="table-wrap"><table>
          <thead><tr>
            <th>Variant</th><th>Genomic</th><th>Consequence</th>
            <th class="num">pLDDT</th><th class="num">RNA LFC</th>
            <th class="num">Splicing</th><th class="num">Access.</th>
          </tr></thead>
          <tbody id="rows">{rows}</tbody>
        </table></div>
        <p class="note">{score_note}</p>
      </div>
    </section>
    <aside>
      <div class="panel">
        <h2>Selected variant</h2>
        <dl id="detail"></dl>
      </div>
      <div class="panel">
        <h2>AlphaFold model</h2>
        <dl>
          <dt>Entry</dt><dd><a href="https://alphafold.ebi.ac.uk/entry/{accession}"
            target="_blank" rel="noopener">{entry_id}</a></dd>
          <dt>Protein</dt><dd>{description}</dd>
          <dt>Mean pLDDT</dt><dd>{mean_plddt}</dd>
        </dl>
      </div>
      <div class="panel legend">
        <h2>Model confidence (pLDDT)</h2>
        {legend}
        <p class="note">pLDDT is confidence in the predicted structure, not the effect of the
        mutation. Low-confidence regions can still be functional (flexible loops, switches).</p>
      </div>
    </aside>
  </div>
  <div class="panel prov">
    <h2>About this snapshot</h2>
    <p class="warn">Research use only, not for clinical decision-making. ClinVar, Ensembl,
    AlphaFold DB and AlphaGenome are updated over time: these results reflect the versions
    below, on the date shown.</p>
    <dl>{provenance}</dl>
    <ul>{attribution}</ul>
  </div>
</main>
<script src="{threedmol_url}"></script>
<script>
(function () {{
  const pdb = {pdb_json};
  const variants = {variants_json};
  const bands = {bands_json};
  const css = getComputedStyle(document.documentElement);
  const accent = css.getPropertyValue("--accent").trim() || "#c2185b";
  const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{
    backgroundColor: css.getPropertyValue("--viewer-bg").trim() || "white",
  }});
  viewer.addModel(pdb, "pdb");
  const plddtColor = (atom) => {{
    for (const [cutoff, colour] of bands) if (atom.b > cutoff) return colour;
    return bands[bands.length - 1][1];
  }};
  const onModel = variants.filter((v) => v.residue !== null);
  const residues = [...new Set(onModel.map((v) => v.residue))];
  let labelsAll = false, spinning = false, surface = null, current = null;

  const fmt = (x) => (x === null ? "-" : (x >= 0 ? "+" : "") + x.toFixed(3));
  const esc = (s) => String(s ?? "-").replace(/[&<>"]/g, (c) =>
    ({{ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }})[c]);

  function draw() {{
    viewer.setStyle({{}}, {{ cartoon: {{ colorfunc: plddtColor }} }});
    viewer.removeAllLabels();
    if (residues.length) {{
      const stick = {{ colorscheme: "magentaCarbon", radius: 0.25 }};
      viewer.addStyle({{ resi: residues }}, {{ stick }});
    }}
    const labelled = labelsAll ? onModel : (current && current.residue !== null ? [current] : []);
    const seen = new Set();
    for (const v of labelled) {{
      if (seen.has(v.residue)) continue;
      seen.add(v.residue);
      viewer.addLabel(v.label, {{
        font: "sans-serif", fontSize: 12, fontColor: "white", backgroundColor: accent,
        backgroundOpacity: 0.9, inFront: true,
      }}, {{ resi: v.residue, atom: "CA" }});
    }}
    viewer.render();
  }}

  function select(i) {{
    current = variants[i];
    document.querySelectorAll("#rows tr").forEach((tr, j) => tr.classList.toggle("sel", j === i));
    const v = current;
    const rows = [
      ["Variant", v.label], ["Genomic", v.variant], ["Consequence", v.consequence],
      ["ClinVar", v.clinvar], ["gnomAD AF", v.gnomad === null ? "-" : v.gnomad.toPrecision(3)],
      ["Residue", v.residue === null ? "not on this model" : v.residue],
      ["pLDDT", v.plddt === null ? "-" : v.plddt.toFixed(1) + " (" + v.band + ")"],
      ["RNA LFC", fmt(v.rna_lfc)], ["Splicing", fmt(v.splicing)],
      ["Accessibility", fmt(v.accessibility)],
    ];
    document.getElementById("detail").innerHTML = rows
      .map(([k, val]) => "<dt>" + esc(k) + "</dt><dd>" + esc(val) + "</dd>").join("");
    draw();
    if (v.residue !== null) {{ viewer.zoomTo({{ resi: v.residue }}, 600); viewer.zoom(0.3, 600); }}
  }}

  document.querySelectorAll("#rows tr").forEach((tr, i) => {{ tr.onclick = () => select(i); }});
  const toggle = (id, fn) => {{
    const el = document.getElementById(id);
    el.onclick = () => el.setAttribute("aria-pressed", String(fn()));
  }};
  document.getElementById("btn-all").onclick = () => viewer.zoomTo({{}}, 600);
  toggle("btn-labels", () => {{ labelsAll = !labelsAll; draw(); return labelsAll; }});
  toggle("btn-spin", () => {{
    spinning = !spinning; viewer.spin(spinning ? "y" : false); return spinning;
  }});
  toggle("btn-surface", () => {{
    if (surface === null) {{
      const style = {{ opacity: 0.5, colorfunc: plddtColor }};
      surface = viewer.addSurface($3Dmol.SurfaceType.VDW, style);
    }} else {{ viewer.removeSurface(surface); surface = null; }}
    viewer.render();
    return surface !== null;
  }});
  window.addEventListener("resize", () => viewer.resize());

  viewer.zoomTo();
  const first = variants.findIndex((v) => v.residue !== null);
  if (first >= 0) select(first); else draw();
}})();
</script>
</body>
</html>
"""


def _readable_consequence(value: str | None) -> str:
    if not value:
        return "-"
    return ", ".join(t.removesuffix("_variant").replace("_", " ") for t in value.split("&"))


def variant_rows(results: Sequence[PipelineResult], structure: StructureResult) -> list[dict]:
    """Per-variant data for the page (first condition's scores; same AlphaFold entry only)."""
    rows: list[dict[str, Any]] = []
    for r in results:
        a = r.annotation
        s = r.structures[0] if r.structures else None
        on_model = s is not None and s.entry.entry_id == structure.entry.entry_id and s.residue
        vs = r.variant_scores[0] if r.variant_scores else None
        rows.append(
            {
                "label": r.vcf_id or r.query,
                "variant": r.query,
                "consequence": _readable_consequence(
                    "&".join(a.consequence_terms)
                    if a and a.consequence_terms
                    else (a.most_severe_consequence if a else None)
                ),
                "clinvar": ", ".join(x.replace("_", " ") for x in a.clinvar_significance)
                if a and a.clinvar_significance
                else None,
                "gnomad": a.gnomad_af if a else None,
                "residue": s.residue if on_model else None,
                "plddt": s.residue_plddt if on_model else None,
                "band": s.residue_confidence if on_model else None,
                "rna_lfc": _top_effect(vs, {"RNA_SEQ"}, "GeneMaskLFC") if vs else None,
                "splicing": _top_effect(
                    vs, {"SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"}
                )
                if vs
                else None,
                "accessibility": _top_effect(vs, {"DNASE", "ATAC"}) if vs else None,
                "condition": vs.label if vs else None,
            }
        )
    return rows


def _num(x: float | None, digits: int = 3, signed: bool = True) -> str:
    if x is None:
        return "-"
    return f"{x:+.{digits}f}" if signed else f"{x:.1f}"


def _provenance_rows(prov: dict[str, Any]) -> str:
    sw, ds = prov.get("software", {}), prov.get("data_sources", {})
    ag = ds.get("alphagenome") or {}
    af = ds.get("alphafold_entries") or {}
    items = [
        ("Generated", prov.get("generated_utc")),
        (
            "geneme-structure",
            f"{sw.get('geneme_structure')} (commit {sw.get('git_commit') or '?'})",
        ),
        ("Assembly", ds.get("genome_assembly")),
        ("Ensembl release", ds.get("ensembl_release")),
        ("ClinVar queried", ds.get("clinvar_query_date")),
        (
            "AlphaFold DB",
            "; ".join(f"{k} v{v.get('model_version') or '?'}" for k, v in af.items()) or None,
        ),
        ("AlphaGenome", ag.get("model")),
        ("AlphaGenome SDK", sw.get("alphagenome_sdk")),
        (
            "Context window",
            ", ".join(f"{w:,} bp" for w in ag.get("sequence_lengths") or []) or None,
        ),
        ("Condition", ", ".join(ag.get("conditions") or []) or None),
    ]
    return "".join(f"<dt>{_e(k)}</dt><dd>{_e(v)}</dd>" for k, v in items if v)


def render_gene_html(
    results: Sequence[PipelineResult],
    *,
    title: str,
    provenance: dict[str, Any],
    nav_html: str = "",
) -> str:
    with_structure = [r for r in results if r.structures]
    if not with_structure:
        raise ValueError("No AlphaFold structure among these results")
    structure = with_structure[0].structures[0]
    if Path(structure.structure_path).suffix != ".pdb":
        raise ValueError("The viewer embeds PDB files; fetch structures with fmt='pdb'")
    entry = structure.entry
    rows = variant_rows(results, structure)

    table = ""
    for v in rows:
        cls = "" if v["residue"] is not None else ' class="off"'
        table += (
            f'<tr{cls}><td><button type="button">{_e(v["label"])}</button></td>'
            f"<td>{_e(v['variant'])}</td><td>{_e(v['consequence'])}</td>"
            f'<td class="num">{_num(v["plddt"], signed=False)}</td>'
            f'<td class="num">{_num(v["rna_lfc"])}</td>'
            f'<td class="num">{_num(v["splicing"])}</td>'
            f'<td class="num">{_num(v["accessibility"])}</td></tr>'
        )
    conditions = {v["condition"] for v in rows if v["condition"]}
    score_note = (
        "AlphaGenome effect scores"
        + (f" in {', '.join(sorted(conditions))}" if conditions else "")
        + ": strongest signed effect per scorer family. Raw scores, no calibrated thresholds:"
        " compare variants down a column, never across columns."
        if conditions
        else "AlphaGenome was not run for these variants."
    )
    legend = "".join(
        f'<div><span style="background:{c}"></span>{_e(t)}</div>' for _, c, t in PLDDT_BANDS
    )
    attribution = "".join(
        f"<li><strong>{_e(k)}.</strong> {_e(v)}</li>"
        for k, v in (provenance.get("attribution") or {}).items()
    )
    n_on = sum(v["residue"] is not None for v in rows)
    gene = entry.gene or (
        with_structure[0].annotation.gene_symbol if with_structure[0].annotation else ""
    )
    return _PAGE.format(
        title=_e(title),
        nav=nav_html,
        heading=_e(title),
        subheading=_e(
            f"{len(rows)} variant(s) in {gene}, {n_on} placed on the AlphaFold model of "
            f"{entry.description or entry.uniprot_accession}. Click a variant to fly to it."
        ),
        entry_id=_e(entry.entry_id),
        accession=_e(entry.uniprot_accession),
        description=_e(entry.description),
        mean_plddt=_e(None if structure.mean_plddt is None else f"{structure.mean_plddt:.1f}"),
        rows=table,
        score_note=_e(score_note),
        legend=legend,
        provenance=_provenance_rows(provenance),
        attribution=attribution,
        threedmol_url=THREEDMOL_URL,
        pdb_json=_script_json(Path(structure.structure_path).read_text()),
        variants_json=_script_json(rows),
        bands_json=_script_json([[cut, colour] for cut, colour, _ in PLDDT_BANDS]),
    )


def write_gene_html(results: Sequence[PipelineResult], path: Path, **kwargs: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_gene_html(results, **kwargs))
    return path
