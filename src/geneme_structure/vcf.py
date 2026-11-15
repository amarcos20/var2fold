"""Minimal, dependency-free VCF reader (plain or gzip/bgzip).

Only the fields needed for scoring are parsed (CHROM, POS, ID, REF, ALT).
Multi-allelic records are split into one variant per ALT. Records that cannot
be scored are skipped with a logged reason: symbolic or breakend ALTs, '*',
non-ACGTN bases, unsupported contigs.
"""

from __future__ import annotations

import gzip
import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from pydantic import ValidationError

from geneme_structure.models import GenomicVariant

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class VcfRecord:
    variant: GenomicVariant
    vcf_id: str | None
    line_no: int


def _open(path: Path):
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def read_vcf(path: str | Path, max_variants: int | None = None) -> Iterator[VcfRecord]:
    """Yield scorable variants from a VCF, in file order."""
    path = Path(path)
    seen = skipped = 0
    with _open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_no}: expected >=5 tab-separated columns")
            chrom, pos, vid, ref, alts = fields[:5]
            for alt in alts.split(","):
                if alt in {".", "*"} or alt.startswith("<") or "[" in alt or "]" in alt:
                    logger.info("%s:%d: skipping non-scorable ALT %r", path.name, line_no, alt)
                    skipped += 1
                    continue
                try:
                    variant = GenomicVariant(
                        chromosome=chrom,
                        position=int(pos),
                        reference_bases=ref,
                        alternate_bases=alt,
                    )
                except (ValidationError, ValueError) as exc:
                    reason = exc.errors()[0]["msg"] if isinstance(exc, ValidationError) else exc
                    logger.warning(
                        "%s:%d: skipping %s:%s %s>%s (%s)",
                        path.name,
                        line_no,
                        chrom,
                        pos,
                        ref,
                        alt,
                        reason,
                    )
                    skipped += 1
                    continue
                yield VcfRecord(
                    variant=variant, vcf_id=None if vid == "." else vid, line_no=line_no
                )
                seen += 1
                if max_variants is not None and seen >= max_variants:
                    logger.warning("Stopped after %d variants (--max-variants)", max_variants)
                    return
    if skipped:
        logger.info("Skipped %d non-scorable allele(s) in %s", skipped, path.name)
