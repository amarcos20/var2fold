"""Write notebooks/gene_report.ipynb, a papermill template that works for any gene.

    python notebooks/build_gene_report.py            # writes the template (not executed)
    papermill notebooks/gene_report.ipynb results/BRCA1_report.ipynb -p gene BRCA1

Run from the repo root.
"""

from pathlib import Path

import nbformat as nbf

md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell

parameters = code(
    """gene = "TP53"                   # HGNC symbol
consequence = "missense"        # missense | coding | splice | any
significance = "pathogenic,likely_pathogenic"
max_variants = 10               # most trusted ClinVar variants; each costs AlphaGenome calls
conditions = ["UBERON:0001157"] # tissue / cell type (CURIEs); first one is used for figures
sequence_length = 1_048_576     # AlphaGenome context window
use_alphagenome = True          # False: annotation + structure only, no key needed
out_dir = "results/gene_reports"
"""
)
parameters.metadata["tags"] = ["parameters"]

cells = [
    md(
        """# Gene report: ClinVar variants → regulation and structure

This notebook runs `genome-structure` for **any human gene**:

1. **ClinVar** gives the known pathogenic variants of the gene (most trusted review status first).
2. **Ensembl VEP** annotates each one (consequence, residue, gnomAD).
3. **AlphaGenome** scores the predicted effect on expression, splicing and accessibility.
4. **AlphaFold DB** places each variant on the protein structure, with **UniProt** domains.

Run it for another gene with [papermill](https://papermill.readthedocs.io):

```bash
papermill notebooks/gene_report.ipynb results/BRCA1_report.ipynb -p gene BRCA1 -p max_variants 15
```

Every number and sentence below is computed from the data. The *How to read this* section at the end explains what the scores do and do not mean."""
    ),
    parameters,
    code(
        """import logging
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import Markdown, display

from genome_structure.alphafold_client import AlphaFoldClient
from genome_structure.clinvar import write_vcf
from genome_structure.config import Settings
from genome_structure.pipeline import GenomeToStructurePipeline, batch_findings, batch_table
from genome_structure.plots import effects_panels, plddt_profile

logging.basicConfig(level=logging.WARNING)
warnings.filterwarnings("ignore", message="IProgress not found")
pd.set_option("display.width", 160)
pd.set_option("display.precision", 3)

ROOT = Path.cwd() if (Path.cwd() / "pyproject.toml").exists() else Path.cwd().parent
os.chdir(ROOT)  # .env, .cache/ and data/ are resolved from the project root
OUT = ROOT / out_dir
OUT.mkdir(parents=True, exist_ok=True)

pipeline = GenomeToStructurePipeline(
    Settings.from_env(), use_alphagenome=use_alphagenome, sequence_length=sequence_length
)"""
    ),
    md("## 1. Known pathogenic variants (ClinVar)"),
    code(
        """selection, results = pipeline.run_gene(
    gene,
    conditions,
    consequence=consequence,
    significance=[s.strip() for s in significance.split(",")],
    max_variants=max_variants,
    score=use_alphagenome,
)
vcf_path = write_vcf(selection, OUT / f"{gene}_clinvar_{consequence}.vcf")
shown = vcf_path.relative_to(ROOT) if vcf_path.is_relative_to(ROOT) else vcf_path
print(f"Selection saved to {shown} (re-run it with: genome-structure vcf ...)")
pd.DataFrame(
    {
        "label": cv.label,
        "variant": str(cv.variant),
        "classification": cv.classification,
        "review stars": cv.stars,
        "consequence": ", ".join(cv.consequences),
        "ClinVar": cv.accession,
    }
    for cv in selection
)"""
    ),
    md("## 2. Annotation, structure and AlphaGenome scores"),
    code(
        """table = batch_table(results)
table.to_csv(OUT / f"{gene}_clinvar_{consequence}_summary.csv", index=False)
cols = ["vcf_id", "variant", "consequence", "protein_change", "gnomad_af", "residue_plddt",
        "rna_lfc", "splicing", "accessibility", "warnings"]
table[[c for c in cols if c in table and table[c].notna().any()]]"""
    ),
    md(
        """## 3. Where do the variants sit on the protein?

AlphaFold's per-residue confidence (pLDDT) along the sequence, the variants analysed, and the UniProt domain architecture underneath."""
    ),
    code(
        """with_structure = [r for r in results if r.structures]
if not with_structure:
    display(Markdown("_No AlphaFold structure for these variants "
                     "(non-coding, or no model in AlphaFold DB)._"))
else:
    s0 = with_structure[0].structures[0]
    plddt = AlphaFoldClient.residue_plddt(Path(s0.structure_path))
    offset = (s0.entry.uniprot_start or 1) - 1  # fragment entries start later in the sequence
    plddt = {res + offset: v for res, v in plddt.items()}
    variants = {
        r.vcf_id or r.query: r.structures[0].residue
        for r in with_structure
        if r.structures[0].entry.entry_id == s0.entry.entry_id and r.structures[0].residue
    }
    domains = pipeline.mapper.protein_features(s0.entry.uniprot_accession)
    fig = plddt_profile(
        plddt, variants, domains,
        title=f"{gene} ({s0.entry.uniprot_accession}) AlphaFold confidence, "
        "with the variants analysed",
    )
    fig.savefig(OUT / f"{gene}_plddt_profile.png", bbox_inches="tight")
    plt.show()"""
    ),
    md(
        """## 4. Predicted regulatory effects (AlphaGenome)

One representative scorer per panel (strongest signed effect over the tracks of the first condition). Scorers have different scales: compare variants **within** a panel, never across panels."""
    ),
    code(
        """scored = table.dropna(subset=["condition"]) if "condition" in table else table.iloc[0:0]
if scored.empty:
    display(Markdown("_AlphaGenome was not run (`use_alphagenome = False`)._"))
else:
    first = scored[scored["condition"] == scored["condition"].iloc[0]]
    window = "1 Mb" if sequence_length == 1_048_576 else f"{sequence_length // 1024} kb"
    fig = effects_panels(
        first.assign(label=first["vcf_id"].fillna(first["variant"])),
        title=f"{gene}: AlphaGenome effect scores, {first['condition'].iloc[0]}, {window} context",
    )
    fig.savefig(OUT / f"{gene}_effects.png", bbox_inches="tight")
    plt.show()"""
    ),
    md("## 5. What the data say"),
    code("""display(Markdown("\\n".join(f"* {line}" for line in batch_findings(table))))"""),
    md(
        """## How to read this

* **pLDDT is confidence in the predicted structure, not the mutation's effect.** A variant in a high-pLDDT region sits in a well-folded part of the protein, where a substitution can disrupt folding or a binding surface. A **low pLDDT does not mean unimportant**: flexible loops can be functional switches (BRAF V600E sits in a low-confidence activation loop).
* **AlphaGenome scores are raw model predictions without calibrated thresholds.** Compare variants within one scorer, not across scorers. Coding missense variants usually get small regulatory scores, because the model sees DNA, not protein function. Splice-region and non-coding variants (`consequence = "splice"` or `"any"`) are where AlphaGenome is most informative.
* **The selection is ranked by ClinVar review status**, then by position, so it is not a random sample of the gene's variants. Raise `max_variants` to cover more of the protein.
* **External databases change** (ClinVar monthly, Ensembl quarterly, AlphaFold DB model versions). Re-running later can give a different selection; the saved VCF lets you re-run the exact same variants.
* This is a research tool for a first look at a set of variants, **not a clinical classifier**."""
    ),
]

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = Path("notebooks/gene_report.ipynb")
nbf.write(nb, out)
print("wrote", out)
