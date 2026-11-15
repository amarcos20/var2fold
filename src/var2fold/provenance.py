"""Record where a set of results came from: date, software and database versions.

External sources change (ClinVar monthly, Ensembl quarterly, AlphaFold DB model
versions, AlphaGenome model updates), so results are only interpretable next to
the versions that produced them. `collect_provenance` gathers what can be
known locally plus the live Ensembl release; every lookup is best-effort.
"""

from __future__ import annotations

import importlib.metadata
import logging
import platform
import subprocess
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from var2fold import __version__
from var2fold.models import PipelineResult

logger = logging.getLogger(__name__)

ENSEMBL_INFO = "https://rest.ensembl.org/info/software"

ATTRIBUTION = {
    "AlphaGenome": (
        "Predictions from AlphaGenome (Google DeepMind), API used under the AlphaGenome "
        "Terms of Service (https://alphagenome.google/terms): outputs are for non-commercial "
        "research use only, not for clinical decision-making, and not for training other "
        "models. Avsec et al., Nature 649, 1206-1218 (2026)."
    ),
    "AlphaFold DB": (
        "Structures from the AlphaFold Protein Structure Database (EMBL-EBI), CC-BY-4.0. "
        "Jumper et al., Nature 2021; Varadi et al., Nucleic Acids Res. 2024."
    ),
    "ClinVar": "Variant classifications from ClinVar (NCBI). Landrum et al., Nucleic Acids Res.",
    "Ensembl": "Annotation from Ensembl VEP and REST. Martin et al., Nucleic Acids Res. 2023.",
    "UniProt": "Protein domains from UniProtKB, CC-BY-4.0.",
}


def _package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def git_commit(repo: Path | None = None) -> str | None:
    """Short commit hash, with '+dirty' when there are uncommitted changes."""
    cwd = repo or Path(__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return f"{commit}+dirty" if dirty else commit


def ensembl_release(session: requests.Session | None = None, timeout: float = 10) -> int | None:
    try:
        resp = (session or requests).get(
            ENSEMBL_INFO, headers={"Accept": "application/json"}, timeout=timeout
        )
        resp.raise_for_status()
        return int(resp.json()["release"])
    except (requests.RequestException, KeyError, ValueError, TypeError) as exc:
        logger.warning("Could not read the Ensembl release: %s", exc)
        return None


def collect_provenance(
    results: Sequence[PipelineResult] = (),
    *,
    session: requests.Session | None = None,
    parameters: dict[str, Any] | None = None,
    clinvar_date: str | None = None,
    commit: str | None = None,
) -> dict[str, Any]:
    """Versions and parameters behind `results` (JSON-serialisable)."""
    entries = {
        s.entry.entry_id: {
            "model_version": s.entry.latest_version,
            "model_created": s.entry.model_created_date,
            "file": Path(s.structure_path).name,
        }
        for r in results
        for s in r.structures
    }
    windows = sorted({vs.interval.width for r in results for vs in r.variant_scores})
    conditions = sorted({vs.condition for r in results for vs in r.variant_scores})
    return {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "software": {
            "var2fold": __version__,
            "git_commit": commit or git_commit(),
            "python": platform.python_version(),
            "alphagenome_sdk": _package_version("alphagenome"),
        },
        "data_sources": {
            "genome_assembly": "GRCh38",
            "ensembl_release": ensembl_release(session),
            "clinvar_query_date": clinvar_date,
            "alphafold_entries": entries,
            "alphagenome": {
                "model": "AlphaGenome 2 (alphagenome.models.v2)" if windows else None,
                "sequence_lengths": windows,
                "conditions": conditions,
            },
        },
        "parameters": parameters or {},
        "attribution": ATTRIBUTION,
    }
