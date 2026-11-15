from pathlib import Path

import pytest

from var2fold.alphafold_client import AlphaFoldClient
from var2fold.models import (
    AlphaFoldEntry,
    PipelineResult,
    StructureResult,
    VariantAnnotation,
)
from var2fold.report import render_structure_html, write_structure_html

# Three-residue toy model; B-factor column holds pLDDT like AlphaFold models.
PDB = """\
ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 95.00           C
ATOM      2  CA  ARG A   2       3.800   0.000   0.000  1.00 72.50           C
ATOM      3  CA  HIS A   3       7.600   0.000   0.000  1.00 40.00           C
END
"""


@pytest.fixture
def structure(tmp_path: Path) -> StructureResult:
    pdb = tmp_path / "AF-P04637-F1-model_v6.pdb"
    pdb.write_text(PDB)
    entry = AlphaFoldEntry.model_validate(
        {
            "entryId": "AF-P04637-F1",
            "uniprotAccession": "P04637",
            "gene": "TP53",
            "uniprotDescription": "Cellular tumor antigen p53 </script><b>x</b>",
            "organismScientificName": "Homo sapiens",
            "uniprotStart": 1,
            "uniprotEnd": 3,
        }
    )
    return StructureResult(
        entry=entry,
        structure_path=str(pdb),
        mean_plddt=69.17,
        residue=2,
        residue_plddt=72.5,
        residue_confidence="confident",
    )


def test_residue_plddt_reads_bfactors(structure):
    assert AlphaFoldClient.residue_plddt(Path(structure.structure_path)) == {
        1: 95.0,
        2: 72.5,
        3: 40.0,
    }


def test_html_embeds_structure_and_escapes_text(structure):
    result = PipelineResult(
        query="chr17:7675088:C>T",
        annotation=VariantAnnotation(gene_symbol="TP53", protein_position=2, amino_acids="R/H"),
        structures=[structure],
    )
    page = render_structure_html(result, structure)
    assert "const residue = 2;" in page
    assert "Focus residue 2" in page
    assert "R/H at residue 2" in page
    assert "chr17:7675088:C&gt;T" in page
    # user/remote text is escaped, never injected as markup
    assert "<b>x</b>" not in page and "&lt;b&gt;x&lt;/b&gt;" in page
    assert page.count("</script>") == 2  # only our own two script tags


def test_write_requires_pdb(structure, tmp_path):
    cif = structure.model_copy(update={"structure_path": str(tmp_path / "x.cif")})
    with pytest.raises(ValueError):
        write_structure_html(PipelineResult(query="q"), cif, tmp_path / "out.html")
