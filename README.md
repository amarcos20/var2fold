# geneme-structure

[![CI](https://github.com/amarcos20/geneme-structure/actions/workflows/ci.yml/badge.svg)](https://github.com/amarcos20/geneme-structure/actions/workflows/ci.yml)

**From a genetic variant to its predicted regulatory effect to the 3D protein it hits, in one command.**

`geneme-structure` links three data sources:

1. **Ensembl VEP / REST** finds the gene, the consequence, the protein position and the UniProt accession for a GRCh38 variant or region. It also returns the allele-specific ClinVar significance, rsID and gnomAD frequency.
2. **[AlphaGenome](https://deepmind.google.com/science/alphagenome)** (Google DeepMind) scores the variant with the model's recommended variant scorers (expression, chromatin accessibility, histone marks, TF binding, methylation, splicing, 3D contacts) over up to 1 Mb of sequence context, using the AlphaGenome 2 SDK.
3. **[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk)** (EMBL-EBI) supplies the predicted 3D structure, with per-residue pLDDT confidence at the mutated residue.

![3D viewer for TP53 R175H: AlphaFold model coloured by pLDDT, Arg175 highlighted, AlphaGenome effect scores alongside](docs/img/tp53_r175h_viewer.png)

## Architecture

```mermaid
flowchart LR
    CLI["CLI<br/>pipeline.py"] --> P[GenomeToStructurePipeline]
    P --> GM["GeneMapper<br/>gene_mapping.py"]
    P --> AG["AlphaGenomeClient<br/>alphagenome_client.py"]
    P --> AF["AlphaFoldClient<br/>alphafold_client.py"]

    GM -->|VEP / overlap| ENS[(Ensembl REST)]
    GM -->|gene → accession| UP[(UniProt REST)]
    AG -->|gRPC via SDK| AGAPI[(AlphaGenome API)]
    AG <-->|sha256-keyed JSON| CACHE[[.cache/alphagenome]]
    AF -->|/api/prediction| AFDB[(AlphaFold DB)]
    AF -->|PDB/mmCIF| FILES[[data/structures]]
    AF -->|Bio.PDB pLDDT| P
    P --> R["ReportWriter<br/>report.py · 3Dmol.js"]

    P --> OUT["results/*.json · *_scores.csv<br/>*.html 3D viewer"]
```

| Module | Responsibility |
|---|---|
| `config.py` | Reads settings strictly from env vars (`.env` via python-dotenv). The API key is held as a `SecretStr`. |
| `models.py` | Pydantic models: `GenomicVariant`, `GenomicInterval`, `VariantAnnotation`, `TrackSummary`, `AlphaFoldEntry`, `PipelineResult` |
| `alphagenome_client.py` | Wraps the official SDK, adds a content-addressed JSON cache and a throttle, and reduces numpy tracks to compact summaries |
| `alphafold_client.py` | AlphaFold DB metadata, cached structure download, per-residue pLDDT parsing |
| `gene_mapping.py` | Variant/region → gene → UniProt, plus ClinVar, rsID and gnomAD AF (batched Ensembl VEP, Ensembl overlap, UniProt search) |
| `vcf.py` | Dependency-free VCF reader (plain or gzip); splits multi-allelic records and skips symbolic ALTs |
| `http_session.py` | Shared `requests` session that retries 429/5xx and honours `Retry-After` |
| `report.py` | Self-contained HTML report: 3Dmol.js viewer coloured by pLDDT, mutated residue highlighted |
| `pipeline.py` | Orchestration, reporting, and the CLI |

### Design notes

- **AlphaGenome has no public REST endpoint.** Access goes through the official `alphagenome` Python SDK (gRPC, `alphagenome.models.v2`). The SDK is an optional extra, imported lazily, so the structure and annotation features work without it.
- **Variant scores vs raw tracks:** the headline numbers come from `score_variant` with AlphaGenome's recommended scorers. These are the model's own recipes, e.g. a log2 ratio in a 501 bp window around the variant, or gene-level log fold change. The summary keeps the strongest score per (output type, scorer) and separates *effect* scorers (ALT vs REF) from *activity* scorers (max of REF/ALT). `--outputs` additionally summarises raw REF/ALT tracks.
- **Windows match the SDK:** `GenomicInterval.resized` reproduces `genome.Interval.resize` exactly. AlphaGenome is not perfectly shift-invariant, so a 1 bp offset in the window would change the scores.
- **Caching:** each AlphaGenome request is canonicalised (kind, variant or region, window, condition, outputs) and hashed with SHA-256. The reduced result is stored as JSON under `.cache/alphagenome/`. Repeated runs cost zero API quota. Pass `--no-cache` to force a live call.
- **Conditions:** AlphaGenome 2 predicts for one *condition*, built from at most one ontology term per type, e.g. `UBERON:0000178_CL:0001054` (blood, CD14+ monocytes). Repeat `-c/--condition` to compare several. Each condition is its own cached call, and the results are pivoted into one table with a column per biosample.
- **Coordinates:** GRCh38 throughout. Variants are 1-based (VCF). Intervals are parsed from 1-based UCSC strings and stored 0-based half-open (BED/AlphaGenome).
- **Secrets:** the key is read only from `os.getenv("ALPHAGENOME_API_KEY")`. The placeholder value counts as unset, and the key is masked in every `repr`.

## Installation

```bash
git clone <your-fork-url> geneme-structure && cd geneme-structure
python -m venv .venv && source .venv/bin/activate
pip install -e ".[alphagenome,dev]"

cp .env.example .env    # then edit .env and set ALPHAGENOME_API_KEY
```

You can request a free non-commercial AlphaGenome API key at <https://deepmind.google.com/science/alphagenome>.

## Usage

> Quote variant strings: an unquoted `>` in `C>T` is a shell redirect.

```bash
# TP53 R175H hotspot: regulatory effect in transverse colon + AlphaFold model
geneme-structure variant 'chr17:7675088:C>T' --condition UBERON:0001157

# ...also summarise raw REF/ALT RNA-seq and DNase tracks
geneme-structure variant 'chr17:7675088:C>T' --outputs RNA_SEQ DNASE

# Compare conditions side by side: colon, liver, blood CD14+ monocytes
geneme-structure variant 'chr17:7675088:C>T' \
    -c UBERON:0001157 -c UBERON:0002107 -c UBERON:0000178_CL:0001054

# Protein-coding genes in a region, their predicted signal, and structures
geneme-structure region chr17:7668421-7687490 --max-genes 5

# Batch: every variant in a VCF -> one summary table (VEP annotation is batched)
geneme-structure vcf examples/tp53_hotspots.vcf --max-variants 25

# Structure only (no API key required)
geneme-structure structure P04637 --residue 175

# Annotation + structure, skipping AlphaGenome
geneme-structure variant 'chr17:7675088:C>T' --no-alphagenome
```

Each run prints a summary and writes:

- `results/<query>.json`: the full validated `PipelineResult`
- `results/<query>_scores.csv`: every AlphaGenome variant score in tidy form (condition, output type, scorer, track, gene, biosample, raw score)
- `results/<query>_conditions.csv`: effect scores pivoted with one column per condition (when there are several)
- `results/<query>_tracks.csv`: raw-track summaries when `--outputs` is given (REF/ALT sums, log2FC, max |Δ|, top track)
- `results/<query>_<AF entry>.html`: an interactive 3D viewer (open in any browser), with the structure embedded and coloured by pLDDT and the mutated residue highlighted
- `data/structures/AF-<acc>-F1-model_v*.pdb`: the downloaded structure
- `vcf` mode: `results/<vcf>_summary.csv`, with one row per (variant, condition). Columns: gene, protein change, ClinVar, gnomAD AF, residue pLDDT, and representative RNA LFC, splicing and accessibility effects. Add `--reports` for per-variant JSON/HTML.

### Python API

```python
from geneme_structure.config import Settings
from geneme_structure.models import GenomicVariant
from geneme_structure.pipeline import GenomeToStructurePipeline

pipeline = GenomeToStructurePipeline(Settings.from_env())
result = pipeline.run_variant(
    GenomicVariant.parse("chr17:7675088:C>T"),
    conditions=["UBERON:0001157", "UBERON:0000178_CL:0001054"],
)
print(result.structures[0].residue_plddt)
for vs in result.variant_scores:
    print(vs.label, [(s.output_type, s.top_score) for s in vs.effects()][:3])
```

### Example: TP53 hotspots

`examples/tp53_hotspots.vcf` holds seven cancer hotspots and the common P72R polymorphism. Abridged output for transverse colon, 16 KB window:

| variant | change | ClinVar | residue pLDDT | RNA LFC | splicing | accessibility |
|---|---|---|---:|---:|---:|---:|
| chr17:7676154:G>C | P72R | pathogenic; benign | 48.3 | -0.064 | 0.127 | 0.008 |
| chr17:7675088:C>T | R175H | (likely) pathogenic | 96.6 | -0.012 | 0.041 | 0.047 |
| chr17:7674872:T>C | Y220C | (likely) pathogenic | 97.7 | -0.045 | **0.521** | 0.025 |
| chr17:7673802:C>T | R273H | (likely) pathogenic | 98.6 | -0.007 | 0.065 | -0.029 |

The hotspots sit in the confidently predicted DNA-binding domain (pLDDT > 94). P72R lies in the disordered proline-rich region. AlphaGenome's regulatory effects are small for most of these coding changes, with Y220C's splicing score as the outlier.

## Configuration

| Variable | Default | Purpose |
|---|---|---|
| `ALPHAGENOME_API_KEY` | (required for predictions) | AlphaGenome API key |
| `GENEME_CACHE_DIR` | `.cache` | AlphaGenome response cache |
| `GENEME_DATA_DIR` | `data` | Downloaded structures |
| `GENEME_HTTP_TIMEOUT` | `30` | REST timeout (seconds) |

## Development

```bash
pytest --cov=geneme_structure   # 51 offline tests, ~91% coverage, no network or API key
ruff check . && ruff format --check .
```

Every HTTP dependency (Ensembl, UniProt, AlphaFold DB) is replaced by a recording `FakeSession` (`tests/fakes.py`), and the AlphaGenome model by small fakes. The suite therefore covers the end-to-end pipeline and CLI, batching, caching and error paths deterministically. CI (`.github/workflows/ci.yml`) runs lint and tests on Python 3.11 to 3.14.

## Limitations

- Human (GRCh38) only.
- The AlphaGenome 2 tidy output has raw scores only (no quantile normalisation), so magnitudes are comparable within a scorer but not across scorers.
- Raw-track summaries (`--outputs`) are deliberately coarse (window sums and max |Δ|). Use the variant scores for interpretation.
- AlphaFold DB models are monomer predictions. pLDDT measures model confidence, not the effect of the mutation.

## Licence & attribution

Code: MIT. AlphaGenome outputs are subject to the [AlphaGenome Terms of Use](https://deepmind.google.com/science/alphagenome). AlphaFold DB data is CC-BY-4.0 (Jumper et al. 2021; Varadi et al. 2024).
