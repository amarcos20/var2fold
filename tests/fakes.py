"""Offline stand-ins for `requests.Session` and canned API payloads."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import requests


@dataclass
class FakeResponse:
    status_code: int = 200
    payload: Any = None
    content: bytes = b""
    text: str = ""

    def json(self) -> Any:
        return self.payload

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)  # type: ignore[arg-type]

    def iter_content(self, chunk_size: int = 1024):
        for i in range(0, len(self.content), chunk_size):
            yield self.content[i : i + chunk_size]

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *exc: object) -> None:
        return None


Handler = Callable[[str, str, dict[str, Any]], FakeResponse]


@dataclass
class FakeSession:
    """Routes requests to handlers by URL substring and records every call."""

    routes: list[tuple[str, Handler]] = field(default_factory=list)
    calls: list[tuple[str, str, dict[str, Any]]] = field(default_factory=list)

    def add(self, url_part: str, handler: Handler | FakeResponse) -> FakeSession:
        h = handler if callable(handler) else (lambda m, u, kw, r=handler: r)
        self.routes.append((url_part, h))
        return self

    def _dispatch(self, method: str, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append((method, url, kwargs))
        for part, handler in self.routes:
            if part in url:
                return handler(method, url, kwargs)
        raise AssertionError(f"Unexpected {method} {url}")

    def get(self, url: str, **kwargs: Any) -> FakeResponse:
        return self._dispatch("GET", url, **kwargs)

    def post(self, url: str, **kwargs: Any) -> FakeResponse:
        return self._dispatch("POST", url, **kwargs)


def vep_echo(records_for: Callable[[str], dict[str, Any]]) -> Handler:
    """VEP handler that answers each posted input line with `records_for(line)`."""

    def handler(method: str, url: str, kwargs: dict[str, Any]) -> FakeResponse:
        lines = kwargs["json"]["variants"]
        return FakeResponse(payload=[{"input": line, **records_for(line)} for line in lines])

    return handler


TOY_PDB = """\
ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 95.00           C
ATOM      2  CA  ARG A   2       3.800   0.000   0.000  1.00 72.50           C
ATOM      3  CA  HIS A   3       7.600   0.000   0.000  1.00 40.00           C
END
"""

AF_ENTRY = {
    "entryId": "AF-P04637-F1",
    "uniprotAccession": "P04637",
    "gene": "TP53",
    "uniprotDescription": "Cellular tumor antigen p53",
    "organismScientificName": "Homo sapiens",
    "uniprotStart": 1,
    "uniprotEnd": 3,
    "globalMetricValue": 69.2,
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.pdb",
    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.cif",
}


def clinvar_doc(
    uid: str,
    spdi: str,
    *,
    protein: str = "R2H",
    gene: str = "TP53",
    classification: str = "Pathogenic",
    review: str = "reviewed by expert panel",
    chrom: str = "17",
    consequences: tuple[str, ...] = ("missense variant",),
) -> dict[str, Any]:
    """Minimal ClinVar esummary document."""
    return {
        "uid": uid,
        "accession": f"VCV{int(uid):09d}",
        "title": f"NM_000546.6({gene}):c.x (p.{protein})",
        "protein_change": protein,
        "genes": [{"symbol": gene}],
        "molecular_consequence_list": list(consequences),
        "germline_classification": {"description": classification, "review_status": review},
        "variation_set": [
            {
                "canonical_spdi": spdi,
                "variation_loc": [
                    {"assembly_name": "GRCh38", "chr": chrom, "start": "1", "stop": "1"},
                    {"assembly_name": "GRCh37", "chr": chrom, "start": "9", "stop": "9"},
                ],
            }
        ],
    }


def eutils(docs: list[dict[str, Any]]) -> Handler:
    """E-utilities handler: esearch returns every uid, esummary the requested docs."""
    by_uid = {d["uid"]: d for d in docs}

    def handler(method: str, url: str, kwargs: dict[str, Any]) -> FakeResponse:
        if "esearch" in url:
            ids = list(by_uid)
            return FakeResponse(payload={"esearchresult": {"count": str(len(ids)), "idlist": ids}})
        wanted = kwargs["data"]["id"].split(",")
        result: dict[str, Any] = {"uids": wanted}
        result.update({u: by_uid[u] for u in wanted})
        return FakeResponse(payload={"result": result})

    return handler
