"""Client for the AlphaFold Protein Structure Database REST API.

Docs: https://alphafold.ebi.ac.uk/api-docs
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Literal

import requests
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from var2fold.http_session import build_session
from var2fold.models import UNIPROT_RE, AlphaFoldEntry, StructureResult, plddt_band

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
StructureFormat = Literal["pdb", "cif"]


class AlphaFoldError(RuntimeError):
    pass


class AlphaFoldClient:
    def __init__(
        self,
        data_dir: Path,
        *,
        timeout: float = 30.0,
        session: requests.Session | None = None,
        base_url: str = ALPHAFOLD_API,
    ) -> None:
        self.structures_dir = Path(data_dir) / "structures"
        self.structures_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.session = session or build_session()
        self.base_url = base_url.rstrip("/")
        self._predictions: dict[str, list[AlphaFoldEntry]] = {}  # per-session memo

    # ------------------------------------------------------------------ metadata

    def get_predictions(self, uniprot_accession: str) -> list[AlphaFoldEntry]:
        """Return all AlphaFold DB models for a UniProt accession ([] if none)."""
        acc = uniprot_accession.strip().upper()
        if not UNIPROT_RE.match(acc):
            raise ValueError(f"{uniprot_accession!r} is not a valid UniProt accession")
        if acc in self._predictions:
            return self._predictions[acc]

        url = f"{self.base_url}/prediction/{acc}"
        logger.info("Fetching AlphaFold DB metadata for %s", acc)
        resp = self.session.get(url, timeout=self.timeout)
        if resp.status_code in (400, 404):
            logger.warning("No AlphaFold DB entry for %s (HTTP %s)", acc, resp.status_code)
            self._predictions[acc] = []
            return []
        resp.raise_for_status()
        data = resp.json()
        if not isinstance(data, list):
            raise AlphaFoldError(f"Unexpected response shape from {url}: {type(data).__name__}")
        entries = [AlphaFoldEntry.model_validate(item) for item in data]
        self._predictions[acc] = entries
        return entries

    def get_entry(
        self, uniprot_accession: str, residue: int | None = None
    ) -> AlphaFoldEntry | None:
        """Pick the best entry, preferring one whose fragment covers `residue`."""
        entries = self.get_predictions(uniprot_accession)
        if not entries:
            return None
        if residue is not None:
            covering = [e for e in entries if e.covers(residue)]
            if covering:
                return covering[0]
            logger.warning("No AlphaFold fragment covers residue %d; using first entry", residue)
        return entries[0]

    # ------------------------------------------------------------------ files

    def download_structure(
        self, entry: AlphaFoldEntry, fmt: StructureFormat = "pdb", *, overwrite: bool = False
    ) -> Path:
        url = entry.pdb_url if fmt == "pdb" else entry.cif_url
        if not url:
            raise AlphaFoldError(f"Entry {entry.entry_id} has no {fmt.upper()} URL")

        dest = self.structures_dir / Path(url).name
        if dest.exists() and dest.stat().st_size > 0 and not overwrite:
            logger.info("Using cached structure %s", dest)
            return dest

        logger.info("Downloading %s", url)
        tmp = dest.with_suffix(dest.suffix + ".part")
        with self.session.get(url, timeout=self.timeout, stream=True) as resp:
            resp.raise_for_status()
            with tmp.open("wb") as fh:
                for chunk in resp.iter_content(chunk_size=64 * 1024):
                    fh.write(chunk)
        tmp.replace(dest)
        return dest

    # ------------------------------------------------------------------ analysis

    @staticmethod
    def residue_plddt(structure_path: Path) -> dict[int, float]:
        """Per-residue pLDDT (AlphaFold stores it in the B-factor column), keyed by residue no."""
        path = Path(structure_path)
        parser = MMCIFParser(QUIET=True) if path.suffix == ".cif" else PDBParser(QUIET=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure(path.stem, str(path))
        scores: dict[int, float] = {}
        model = next(structure.get_models())
        for residue in model.get_residues():
            if "CA" in residue:
                scores[residue.id[1]] = float(residue["CA"].get_bfactor())
        return scores

    def fetch_structure(
        self,
        uniprot_accession: str,
        residue: int | None = None,
        fmt: StructureFormat = "pdb",
    ) -> StructureResult | None:
        """Metadata + download + pLDDT summary (optionally at a specific residue)."""
        entry = self.get_entry(uniprot_accession, residue)
        if entry is None:
            return None
        path = self.download_structure(entry, fmt)
        plddt = self.residue_plddt(path)
        mean = sum(plddt.values()) / len(plddt) if plddt else None

        res_score = None
        if residue is not None and entry.uniprot_start is not None:
            # Model residue numbering starts at 1 for the fragment's first residue.
            res_score = plddt.get(residue - entry.uniprot_start + 1)
        elif residue is not None:
            res_score = plddt.get(residue)

        return StructureResult(
            entry=entry,
            structure_path=str(path),
            mean_plddt=round(mean, 2) if mean is not None else None,
            residue=residue,
            residue_plddt=res_score,
            residue_confidence=plddt_band(res_score) if res_score is not None else None,
        )
