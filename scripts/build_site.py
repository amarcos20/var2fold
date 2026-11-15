"""Build the GitHub Pages site in docs/ (run from the repo root).

    python scripts/build_site.py

For each case study it writes docs/<slug>/viewer.html (all variants on the
AlphaFold model), docs/<slug>/report.html (the executed notebook) and
docs/<slug>/provenance.json, plus docs/index.html. Variants come from fixed,
dated VCFs in examples/, and AlphaGenome responses come from the local cache,
so rebuilding is reproducible and costs no API quota once the examples have run.
"""

from __future__ import annotations

import html
import json
import logging
import re
from datetime import date
from pathlib import Path

import nbformat
from nbconvert import HTMLExporter

from var2fold.config import Settings
from var2fold.gene_viewer import write_gene_html
from var2fold.pipeline import GenomeToStructurePipeline
from var2fold.provenance import ATTRIBUTION, collect_provenance, git_commit

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
COLON = "UBERON:0001157"
WINDOW = 1_048_576
REPO_URL = "https://github.com/amarcos20/var2fold"

CASES = [
    {
        "slug": "brca1",
        "title": "BRCA1: pathogenic missense variants from ClinVar",
        "vcf": "examples/BRCA1_clinvar_missense.vcf",
        "notebook": "notebooks/examples/BRCA1_report.ipynb",
        "image": "img/brca1_effects.png",
        "command": "var2fold gene BRCA1 --max-variants 10",
        "blurb": (
            "Ten expert-panel pathogenic missense variants, most in the folded RING and BRCT "
            "domains. Two of them, D1692H and E1559K, get AlphaGenome splicing scores about "
            "40× the median. Both sit on the last base of an exon, and VEP independently flags "
            "exactly these two as splice-region variants: damage that may come from "
            "mis-splicing as well as from the amino-acid change."
        ),
    },
    {
        "slug": "tp53",
        "title": "TP53: cancer hotspots",
        "vcf": "examples/tp53_hotspots.vcf",
        "notebook": "notebooks/tp53_hotspots.ipynb",
        "image": "img/tp53_plddt_profile.png",
        "command": "var2fold vcf examples/tp53_hotspots.vcf",
        "blurb": (
            "Seven cancer hotspots and the common P72R polymorphism. Every hotspot sits in the "
            "folded DNA-binding domain (residue pLDDT 95–99), which fits a structural "
            "mechanism. Predicted regulatory effects are small, with one splicing outlier "
            "(Y220C) worth testing."
        ),
    },
]

BENCHMARK = {
    "slug": "benchmark",
    "title": "Benchmark: 300 ClinVar missense variants",
    "notebook": "notebooks/clinvar_benchmark.ipynb",
    "image": "img/benchmark_roc.png",
    "command": "python notebooks/build_benchmark.py",
    "blurb": (
        "Do the scores carry information? On 150 pathogenic and 150 benign expert-reviewed "
        "missense variants in six genes, AlphaFold pLDDT at the mutated residue separates the "
        "classes with AUC 0.82. AlphaGenome splicing is weaker overall (AUC 0.67) but its top "
        "20 scores are 19 pathogenic, and it flags 10 of the 31 pathogenic variants that "
        "structure cannot explain. Expression and accessibility show no signal, as expected "
        "for missense variants."
    ),
}

NAV = (
    '<nav><a href="../">← var2fold</a> · <a href="report.html">Full report</a>'
    ' · <a href="provenance.json">provenance.json</a></nav>'
)


def vcf_file_date(path: Path) -> str | None:
    for line in path.read_text().splitlines():
        if m := re.match(r"##fileDate=(\d{4})(\d{2})(\d{2})", line):
            return "-".join(m.groups())
        if not line.startswith("##"):
            break
    return None


def notebook_html(path: Path, title: str, viewer: bool = True) -> str:
    nb = nbformat.read(path, as_version=4)
    exporter = HTMLExporter(template_name="lab")
    body, _ = exporter.from_notebook_node(nb, resources={"metadata": {"name": title}})
    banner = (
        '<div style="font:14px system-ui,sans-serif;padding:10px 16px;border-bottom:1px solid '
        '#ccc"><a href="../">← var2fold</a> · '
        + ('<a href="viewer.html">3D viewer</a> · ' if viewer else "")
        + "Executed notebook, rendered as a static page. Research use only.</div>"
    )
    return re.sub(r"(<body[^>]*>)", lambda m: m.group(1) + banner, body, count=1)


INDEX = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>var2fold</title>
<meta name="description" content="From a genetic variant to its predicted regulatory effect
and the 3D protein it hits: AlphaGenome + AlphaFold DB + ClinVar.">
<style>
  :root {{
    --bg: #ffffff; --fg: #1d2330; --muted: #5d6676; --panel: #f4f6f9;
    --border: #dde2ea; --accent: #c2185b;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #12151b; --fg: #e6e9ef; --muted: #9aa3b2; --panel: #1b2029;
      --border: #2c3340; --accent: #f06292;
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: var(--bg); color: var(--fg);
    font: 16px/1.6 system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  main {{ max-width: 980px; margin: 0 auto; padding: 40px 16px 56px; }}
  h1 {{ font-size: 2rem; margin: 0 0 6px; letter-spacing: -.01em; }}
  .lead {{ font-size: 1.15rem; color: var(--muted); margin: 0 0 28px; max-width: 60ch; }}
  h2 {{ font-size: 1.15rem; margin: 40px 0 12px; }}
  .flow {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; }}
  @media (max-width: 760px) {{ .flow {{ grid-template-columns: 1fr 1fr; }} }}
  .step {{ background: var(--panel); border: 1px solid var(--border); border-radius: 10px;
    padding: 12px 14px; font-size: .9rem; }}
  .step b {{ display: block; margin-bottom: 2px; }}
  .step span {{ color: var(--muted); }}
  .case {{ display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.1fr); gap: 20px;
    background: var(--panel); border: 1px solid var(--border); border-radius: 12px;
    padding: 18px; margin-bottom: 18px; }}
  @media (max-width: 760px) {{ .case {{ grid-template-columns: 1fr; }} }}
  .case img {{ width: 100%; height: auto; border-radius: 8px; background: #fff;
    border: 1px solid var(--border); }}
  .case h3 {{ margin: 0 0 8px; font-size: 1.1rem; }}
  .case p {{ margin: 0 0 12px; font-size: .95rem; }}
  .links {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 4px; }}
  .links a {{ text-decoration: none; padding: 6px 12px; border-radius: 6px;
    border: 1px solid var(--border); background: var(--bg); font-size: .9rem; }}
  .links a.primary {{ background: var(--accent); border-color: var(--accent); color: #fff; }}
  code, pre {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .85rem; }}
  pre {{ background: var(--panel); border: 1px solid var(--border); border-radius: 8px;
    padding: 12px 14px; overflow-x: auto; }}
  .case code {{ color: var(--muted); }}
  .notice {{ border-left: 3px solid var(--accent); padding: 2px 0 2px 12px; }}
  ul.attr {{ font-size: .85rem; color: var(--muted); padding-left: 18px; }}
  footer {{ color: var(--muted); font-size: .85rem; margin-top: 40px;
    border-top: 1px solid var(--border); padding-top: 16px; }}
  a {{ color: inherit; }}
</style>
</head>
<body>
<main>
  <h1>var2fold</h1>
  <p class="lead">From a genetic variant to its predicted regulatory effect and the 3D
  protein it hits. Start from a variant, a VCF or just a gene name.</p>

  <div class="flow">
    <div class="step"><b>1 · ClinVar</b><span>known pathogenic variants for a gene</span></div>
    <div class="step"><b>2 · Ensembl VEP</b><span>gene, consequence, residue, gnomAD</span></div>
    <div class="step"><b>3 · AlphaGenome</b><span>predicted effect on expression,
      splicing, chromatin</span></div>
    <div class="step"><b>4 · AlphaFold DB</b><span>3D model and confidence at the
      mutated residue</span></div>
  </div>

  <h2>Case studies</h2>
  {cases}

  <h2>Does it work?</h2>
  {benchmark}

  <h2>Run it on your gene</h2>
<pre>pip install "var2fold[alphagenome] @ git+https://github.com/amarcos20/var2fold"
var2fold gene YOUR_GENE --max-variants 10
papermill notebooks/gene_report.ipynb YOUR_GENE.ipynb -p gene YOUR_GENE</pre>
  <p>Code, documentation and tests: <a href="{repo}">{repo}</a></p>

  <h2>About this snapshot</h2>
  <p class="notice">Research use only, not for clinical decision-making. These pages were
  generated on <strong>{generated}</strong>. ClinVar, Ensembl, AlphaFold DB and AlphaGenome
  change over time; each case lists the exact versions in its <code>provenance.json</code>
  and in the panel at the bottom of its 3D viewer.</p>
  <ul class="attr">{attribution}</ul>
  <footer>var2fold {version} · built with <code>scripts/build_site.py</code></footer>
</main>
</body>
</html>
"""

CASE = """<article class="case">
  <img src="{image}" alt="{alt}">
  <div>
    <h3>{title}</h3>
    <p>{blurb}</p>
    <p><code>{command}</code></p>
    <div class="links">
      <a class="primary" href="{slug}/viewer.html">Open 3D viewer</a>
      <a href="{slug}/report.html">Full report</a>
      <a href="{slug}/provenance.json">Provenance</a>
    </div>
  </div>
</article>"""


BENCHMARK_CARD = """<article class="case">
  <img src="{image}" alt="{alt}">
  <div>
    <h3>{title}</h3>
    <p>{blurb}</p>
    <p><code>{command}</code></p>
    <div class="links">
      <a class="primary" href="{slug}/report.html">Full analysis</a>
    </div>
  </div>
</article>"""


def main() -> None:
    logging.basicConfig(level=logging.WARNING)
    settings = Settings.from_env()
    pipeline = GenomeToStructurePipeline(settings, sequence_length=WINDOW)
    commit = git_commit()  # before writing: the build itself modifies tracked files
    cards = []
    for case in CASES:
        out = DOCS / case["slug"]
        out.mkdir(parents=True, exist_ok=True)
        vcf = ROOT / case["vcf"]
        results = pipeline.run_vcf(vcf, conditions=[COLON], max_variants=None)
        provenance = collect_provenance(
            results,
            session=pipeline.mapper.session,
            clinvar_date=vcf_file_date(vcf),
            commit=commit,
            parameters={"vcf": case["vcf"], "conditions": [COLON], "sequence_length": WINDOW},
        )
        (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        write_gene_html(
            results,
            out / "viewer.html",
            title=case["title"],
            provenance=provenance,
            nav_html=NAV,
        )
        (out / "report.html").write_text(notebook_html(ROOT / case["notebook"], case["title"]))
        cards.append(
            CASE.format(
                slug=case["slug"],
                image=case["image"],
                alt=html.escape(f"Figure from the {case['title']} case study"),
                title=html.escape(case["title"]),
                blurb=html.escape(case["blurb"]),
                command=html.escape(case["command"]),
            )
        )
        print(f"wrote docs/{case['slug']}/ ({len(results)} variants)")

    bench_dir = DOCS / BENCHMARK["slug"]
    bench_dir.mkdir(parents=True, exist_ok=True)
    (bench_dir / "report.html").write_text(
        notebook_html(ROOT / BENCHMARK["notebook"], BENCHMARK["title"], viewer=False)
    )
    benchmark_card = BENCHMARK_CARD.format(
        slug=BENCHMARK["slug"],
        image=BENCHMARK["image"],
        alt="ROC curves of each var2fold feature on the ClinVar missense benchmark",
        title=html.escape(BENCHMARK["title"]),
        blurb=html.escape(BENCHMARK["blurb"]),
        command=html.escape(BENCHMARK["command"]),
    )
    print("wrote docs/benchmark/")

    attribution = "".join(
        f"<li><strong>{html.escape(k)}.</strong> {html.escape(v)}</li>"
        for k, v in ATTRIBUTION.items()
    )
    from var2fold import __version__

    (DOCS / "index.html").write_text(
        INDEX.format(
            cases="\n".join(cards),
            benchmark=benchmark_card,
            repo=REPO_URL,
            generated=date.today().isoformat(),
            attribution=attribution,
            version=__version__,
        )
    )
    (DOCS / ".nojekyll").write_text("")  # serve files as-is (no Jekyll processing)
    print("wrote docs/index.html")


if __name__ == "__main__":
    main()
