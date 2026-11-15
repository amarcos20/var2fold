"""Build and execute notebooks/tp53_hotspots.ipynb (run from the repo root)."""

import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell

cells = [
    md("""# TP53 cancer hotspots: structure or regulation?

TP53 is the most frequently mutated gene in human cancer. Most of its somatic mutations are missense changes clustered at a few *hotspot* residues (R175, Y220, G245, R248, R273, R282).

This notebook asks one question with `geneme-structure`: **do these hotspots act through the protein structure, through gene regulation, or both?**

* **AlphaFold DB** gives the 3D model of p53 and its per-residue confidence (pLDDT). Hotspots in a well-folded domain point to a structural mechanism.
* **AlphaGenome** scores each variant's predicted effect on expression, splicing and chromatin accessibility in the full 1 Mb genomic context.
* **Ensembl VEP / ClinVar / gnomAD** supply the clinical and population context.

The common P72R polymorphism (rs1042522) is included as a contrast.

> Live AlphaGenome calls need `ALPHAGENOME_API_KEY` in `.env`. Responses are cached under `.cache/`, so re-running costs nothing."""),
    code("""import logging
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from geneme_structure.alphafold_client import AlphaFoldClient
from geneme_structure.alphagenome_client import compare_conditions
from geneme_structure.config import Settings
from geneme_structure.models import GenomicVariant
from geneme_structure.pipeline import GenomeToStructurePipeline, batch_table, write_report
from geneme_structure.plots import effects_panels, plddt_profile

logging.basicConfig(level=logging.WARNING)
warnings.filterwarnings("ignore", message="IProgress not found")  # tqdm without ipywidgets
pd.set_option("display.width", 160)
pd.set_option("display.precision", 4)

ROOT = Path.cwd() if (Path.cwd() / "pyproject.toml").exists() else Path.cwd().parent
os.chdir(ROOT)  # .env, .cache/ and data/ are resolved from the project root
VCF = ROOT / "examples" / "tp53_hotspots.vcf"
COLON = "UBERON:0001157"  # transverse colon

pipeline = GenomeToStructurePipeline(Settings.from_env(), sequence_length=1_048_576)"""),
    md("""## 1. Annotate and score every variant

One call runs the whole pipeline for the VCF: batched VEP annotation, AlphaGenome recommended scorers (1 Mb window, transverse colon) and the AlphaFold model with the pLDDT at each mutated residue."""),
    code("""results = pipeline.run_vcf(VCF, conditions=[COLON])
table = batch_table(results)
table["label"] = [r.vcf_id for r in results]
table[["label", "variant", "protein_change", "clinvar", "gnomad_af", "residue_plddt",
       "rna_lfc", "splicing", "accessibility"]]"""),
    md("""## 2. Where do the hotspots sit on the structure?

AlphaFold writes pLDDT into the B-factor column of the model. Plotting it along the sequence shows p53's domain architecture: a disordered N-terminal transactivation and proline-rich region, the folded DNA-binding domain (102–292) and the oligomerisation domain (325–356). The bands are AlphaFold DB's confidence thresholds; the domain bars come from UniProt."""),
    code("""structure = results[1].structures[0]  # R175H; all variants share the p53 model
plddt = AlphaFoldClient.residue_plddt(Path(structure.structure_path))
hotspots = {r.vcf_id: r.structures[0].residue for r in results if r.structures}
domains = pipeline.mapper.protein_features(structure.entry.uniprot_accession)  # UniProt

fig = plddt_profile(plddt, hotspots, domains,
                    title="p53 AlphaFold confidence along the sequence, with the variants analysed")
fig.savefig(ROOT / "docs" / "img" / "tp53_plddt_profile.png", bbox_inches="tight")
plt.show()"""),
    md("""## 3. How large are the predicted regulatory effects?

Each panel shows one representative AlphaGenome scorer (strongest signed effect over the colon tracks). The scorers have different scales, so compare variants **within** a panel, never across panels."""),
    code("""ordered = table.sort_values("variant")  # genomic order
fig = effects_panels(ordered, title="AlphaGenome effect scores, transverse colon, 1 Mb context")
fig.savefig(ROOT / "docs" / "img" / "tp53_effects.png", bbox_inches="tight")
plt.show()"""),
    md("""## 4. Is the hotspot's regulatory effect tissue-specific?

AlphaGenome 2 predicts for one *condition* at a time. Here R175H is compared across transverse colon, liver and blood CD14+ monocytes."""),
    code("""r175h = GenomicVariant.parse("chr17:7675088:C>T")
conditions = [COLON, "UBERON:0002107", "UBERON:0000178_CL:0001054"]
multi = pipeline.run_variant(r175h, conditions=conditions)
compare_conditions(multi.variant_scores)"""),
    md("""## 5. Interactive 3D report

`write_report` produces a self-contained HTML page for each structure. It shows the AlphaFold model coloured by pLDDT, R175 highlighted, and the comparison table alongside. Open it in a browser:"""),
    code("""html = [p for p in write_report(multi, ROOT / "results" / "notebook").parent.glob("*.html")]
print(*[p.relative_to(ROOT) for p in html], sep="\\n")"""),
    md("""![3D report for R175H](../docs/img/tp53_r175h_viewer.png)"""),
    md("""## What the data say

**Structure: every hotspot sits in the folded core.** All seven hotspot residues fall in the DNA-binding domain, where AlphaFold is very confident (residue pLDDT 94.8–98.6). The P72R polymorphism lies in the disordered proline-rich region (pLDDT 48). The data fit the established view that hotspots disrupt p53's folding or DNA contacts rather than its expression.

**Regulation: predicted effects are small.** Across the hotspots, AlphaGenome's gene-level RNA-seq log fold changes stay within ±0.014. Accessibility changes stay within ±0.051. For coding missense variants that is the expected result: the model sees DNA sequence, not protein function.

**One outlier worth following up: Y220C and splicing.** Y220C's strongest splicing score (+0.50) is more than twice that of any other variant (next: P72R and R248Q at +0.23). This is a model prediction, not evidence. It is a concrete hypothesis to check against RNA-seq from Y220C carriers or a minigene assay.

**R175H looks similar across tissues.** Its splice-junction score is consistent in colon, liver and CD14+ monocytes (+0.11 to +0.14). The largest between-tissue spread is in methylation (−0.02 to −0.14).

### Caveats
* AlphaGenome 2 reports raw scores without calibrated thresholds. Compare variants within a scorer, not across scorers.
* Window size matters. A 16 kb context does not cover TP53's full gene body, and it gave different (generally smaller) splicing scores for these variants than the 1 Mb context used here.
* pLDDT is the model's confidence in the predicted structure, not a measure of the mutation's effect."""),
]

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = Path("notebooks/tp53_hotspots.ipynb")
out.parent.mkdir(exist_ok=True)
if "--no-exec" not in sys.argv:
    NotebookClient(
        nb, timeout=1800, kernel_name="python3", resources={"metadata": {"path": "notebooks"}}
    ).execute()
nbf.write(nb, out)
print("wrote", out)
