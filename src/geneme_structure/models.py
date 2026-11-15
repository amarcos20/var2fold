"""Pydantic data models shared across clients and the pipeline.

Coordinate conventions
----------------------
* `GenomicVariant.position` is **1-based** (VCF style).
* `GenomicInterval` is **0-based, half-open** `[start, end)` (BED / AlphaGenome style).
* Assembly is GRCh38/hg38 throughout (AlphaGenome and Ensembl REST default).
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}
_BASES = re.compile(r"^[ACGTN]+$")
_VARIANT_RE = re.compile(
    r"^(?P<chrom>[\w.]+)[:_-](?P<pos>\d+)[:_-](?P<ref>[ACGTNacgtn]+)[>:_/-](?P<alt>[ACGTNacgtn]+)$"
)
_INTERVAL_RE = re.compile(r"^(?P<chrom>[\w.]+):(?P<start>[\d,]+)-(?P<end>[\d,]+)$")
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")

UniProtAccession = Annotated[str, Field(pattern=UNIPROT_RE.pattern)]


def normalize_chromosome(value: str) -> str:
    """Return a UCSC-style chromosome name ('chr17', 'chrX', 'chrM')."""
    core = value.strip().removeprefix("chr").removeprefix("Chr").upper()
    if core == "MT":
        core = "M"
    if core not in _VALID_CHROMS:
        raise ValueError(f"Unsupported chromosome: {value!r}")
    return f"chr{core}"


class _Base(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class GenomicVariant(_Base):
    chromosome: str
    position: int = Field(ge=1, description="1-based position of the first REF base")
    reference_bases: str
    alternate_bases: str

    @field_validator("chromosome")
    @classmethod
    def _normalize_chrom(cls, v: str) -> str:
        return normalize_chromosome(v)

    @field_validator("reference_bases", "alternate_bases")
    @classmethod
    def _check_bases(cls, v: str) -> str:
        v = v.upper()
        if not _BASES.match(v):
            raise ValueError(f"Bases must match [ACGTN]+, got {v!r}")
        return v

    @model_validator(mode="after")
    def _ref_ne_alt(self) -> GenomicVariant:
        if self.reference_bases == self.alternate_bases:
            raise ValueError("REF and ALT must differ")
        return self

    @classmethod
    def parse(cls, text: str) -> GenomicVariant:
        """Parse 'chr17:7675088:C>T' (also accepts '17-7675088-C-T')."""
        m = _VARIANT_RE.match(text.strip())
        if not m:
            raise ValueError(f"Cannot parse variant {text!r}; expected e.g. chr17:7675088:C>T")
        return cls(
            chromosome=m["chrom"],
            position=int(m["pos"]),
            reference_bases=m["ref"],
            alternate_bases=m["alt"],
        )

    @property
    def reference_interval(self) -> GenomicInterval:
        """0-based half-open span of the REF allele."""
        return GenomicInterval(
            chromosome=self.chromosome,
            start=self.position - 1,
            end=self.position - 1 + len(self.reference_bases),
        )

    @property
    def ensembl_chromosome(self) -> str:
        return self.chromosome.removeprefix("chr").replace("M", "MT")

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.reference_bases}>{self.alternate_bases}"


class GenomicInterval(_Base):
    chromosome: str
    start: int = Field(ge=0, description="0-based inclusive start")
    end: int = Field(gt=0, description="0-based exclusive end")

    @field_validator("chromosome")
    @classmethod
    def _normalize_chrom(cls, v: str) -> str:
        return normalize_chromosome(v)

    @model_validator(mode="after")
    def _check_order(self) -> GenomicInterval:
        if self.end <= self.start:
            raise ValueError("end must be greater than start")
        return self

    @classmethod
    def parse(cls, text: str) -> GenomicInterval:
        """Parse a 1-based inclusive UCSC string 'chr17:7,668,421-7,687,490'."""
        m = _INTERVAL_RE.match(text.strip())
        if not m:
            raise ValueError(f"Cannot parse interval {text!r}; expected e.g. chr17:7668421-7687490")
        start_1b = int(m["start"].replace(",", ""))
        end_1b = int(m["end"].replace(",", ""))
        return cls(chromosome=m["chrom"], start=start_1b - 1, end=end_1b)

    @property
    def width(self) -> int:
        return self.end - self.start

    @property
    def center(self) -> int:
        """Centre rounded up for odd widths (same convention as `alphagenome.data.genome`)."""
        return (self.start + self.end + 1) // 2

    def resized(self, width: int) -> GenomicInterval:
        """Resize around the centre exactly like `genome.Interval.resize` (+ strand).

        Matching the SDK matters: AlphaGenome is not perfectly shift-invariant,
        so a 1 bp offset in the window changes the scores slightly.
        """
        start = self.center - (width + 1) // 2
        if start < 0:
            raise ValueError(f"{width} bp window around {self} runs off the chromosome start")
        return GenomicInterval(chromosome=self.chromosome, start=start, end=start + width)

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start + 1}-{self.end}"


# --------------------------------------------------------------------------- annotation


class VariantAnnotation(_Base):
    """Gene/protein context for a variant (from Ensembl VEP)."""

    gene_symbol: str | None = None
    gene_id: str | None = None
    transcript_id: str | None = None
    most_severe_consequence: str | None = None
    consequence_terms: list[str] = Field(default_factory=list)
    protein_position: int | None = None
    amino_acids: str | None = None
    uniprot_accessions: list[str] = Field(default_factory=list)

    @property
    def primary_uniprot(self) -> str | None:
        return self.uniprot_accessions[0] if self.uniprot_accessions else None


# --------------------------------------------------------------------------- AlphaGenome


class TrackSummary(_Base):
    """Reduced, JSON-serialisable summary of one AlphaGenome output type."""

    output_type: str
    n_tracks: int
    resolution_bp: int | None = None
    ref_sum: float | None = None
    alt_sum: float | None = None
    log2_fold_change: float | None = None
    max_abs_delta: float | None = None
    top_track: str | None = Field(default=None, description="Track with the largest |ALT-REF|")
    mean_signal: float | None = Field(default=None, description="Region mode only")


class AlphaGenomePrediction(_Base):
    query: str
    interval: GenomicInterval
    ontology_terms: list[str]
    tracks: list[TrackSummary]
    from_cache: bool = False
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class ScoreRecord(_Base):
    """One row of AlphaGenome's tidy variant-score table (selected columns)."""

    output_type: str
    variant_scorer: str
    track_name: str | None = None
    gene_name: str | None = None
    biosample_name: str | None = None
    raw_score: float


class ScorerSummary(_Base):
    """Strongest signal for one (output type, scorer) pair."""

    output_type: str
    scorer: str
    kind: str = Field(description="'effect' (ALT vs REF) or 'activity' (max of REF/ALT)")
    n_scores: int
    max_abs_score: float
    top_score: float
    top_track: str | None = None
    top_gene: str | None = None


class VariantScores(_Base):
    """Output of `score_variant` with AlphaGenome's recommended scorers."""

    query: str
    interval: GenomicInterval
    condition: str
    summary: list[ScorerSummary]
    records: list[ScoreRecord] = Field(default_factory=list, repr=False)
    from_cache: bool = False
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    def effects(self) -> list[ScorerSummary]:
        return [s for s in self.summary if s.kind == "effect"]


# --------------------------------------------------------------------------- AlphaFold


class AlphaFoldEntry(BaseModel):
    """Subset of an AlphaFold DB `/api/prediction/{accession}` record.

    The upstream schema evolves, so unknown fields are ignored and most known
    ones are optional.
    """

    model_config = ConfigDict(extra="ignore", populate_by_name=True, frozen=True)

    entry_id: str = Field(alias="entryId")
    uniprot_accession: str = Field(alias="uniprotAccession")
    uniprot_id: str | None = Field(default=None, alias="uniprotId")
    gene: str | None = None
    description: str | None = Field(default=None, alias="uniprotDescription")
    organism: str | None = Field(default=None, alias="organismScientificName")
    uniprot_start: int | None = Field(default=None, alias="uniprotStart")
    uniprot_end: int | None = Field(default=None, alias="uniprotEnd")
    sequence: str | None = Field(default=None, alias="uniprotSequence")
    latest_version: int | None = Field(default=None, alias="latestVersion")
    model_created_date: str | None = Field(default=None, alias="modelCreatedDate")
    global_plddt: float | None = Field(default=None, alias="globalMetricValue")
    pdb_url: str | None = Field(default=None, alias="pdbUrl")
    cif_url: str | None = Field(default=None, alias="cifUrl")
    pae_image_url: str | None = Field(default=None, alias="paeImageUrl")

    def covers(self, residue: int) -> bool:
        if self.uniprot_start is None or self.uniprot_end is None:
            return True
        return self.uniprot_start <= residue <= self.uniprot_end


class StructureResult(_Base):
    entry: AlphaFoldEntry
    structure_path: str
    mean_plddt: float | None = None
    residue: int | None = None
    residue_plddt: float | None = None
    residue_confidence: str | None = None


def plddt_band(score: float) -> str:
    """AlphaFold DB confidence bands."""
    if score > 90:
        return "very high"
    if score > 70:
        return "confident"
    if score > 50:
        return "low"
    return "very low"


# --------------------------------------------------------------------------- pipeline


class PipelineResult(_Base):
    query: str
    annotation: VariantAnnotation | None = None
    genes_in_region: list[str] = Field(default_factory=list)
    alphagenome: AlphaGenomePrediction | None = None
    variant_scores: VariantScores | None = None
    structures: list[StructureResult] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
