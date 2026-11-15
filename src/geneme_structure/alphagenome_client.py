"""AlphaGenome client wrapper with on-disk JSON caching.

AlphaGenome is served through Google DeepMind's official Python SDK
(``pip install alphagenome``, AlphaGenome 2 / ``alphagenome.models.v2``), which
talks gRPC rather than a public REST endpoint. This module wraps that SDK and adds:

* a content-addressed JSON cache, so repeated queries never hit the API
  (AlphaGenome has per-key quotas);
* a minimum interval between live calls (simple client-side throttle);
* reduction of the large numpy track arrays into compact `TrackSummary`
  records that are easy to cache, log and tabulate.

The SDK is imported lazily so the rest of the package works without it.

Conditions: AlphaGenome 2 predicts for a single *condition*, a composite of at
most one ontology term per type (UBERON tissue, CL cell type, CLO, EFO, NTR),
serialised as ``'UBERON:0000178_CL:0001054'``. ``ontology_terms`` passed to
this client are joined into that condition.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, TypeVar

import numpy as np
import pandas as pd
from pydantic import BaseModel

from geneme_structure.models import (
    AlphaGenomePrediction,
    GenomicInterval,
    GenomicVariant,
    ScoreRecord,
    ScorerSummary,
    TrackSummary,
    VariantScores,
)

logger = logging.getLogger(__name__)

ModelT = TypeVar("ModelT", bound=BaseModel)

# Context windows accepted by AlphaGenome 2 (16 KB, 100 KB, 500 KB, 1 MB).
SUPPORTED_SEQUENCE_LENGTHS: tuple[int, ...] = (16_384, 131_072, 524_288, 1_048_576)

# Output types whose result is a 1D-per-track `TrackData` (positions x tracks).
TRACK_OUTPUT_TYPES: tuple[str, ...] = (
    "RNA_SEQ",
    "CAGE",
    "DNASE",
    "ATAC",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "METHYLATION",
)

CACHE_SCHEMA_VERSION = 3  # bump when the SDK/model or summary format changes


class AlphaGenomeUnavailableError(ImportError):
    """Raised when the `alphagenome` SDK is not installed."""


def to_condition(ontology_terms: Sequence[str]) -> str:
    """Join ontology CURIEs into an AlphaGenome 2 composite condition string."""
    terms = [t.strip() for t in ontology_terms if t.strip()]
    if not terms:
        raise ValueError("At least one ontology term (e.g. UBERON:0001157) is required")
    prefixes = [t.split(":", 1)[0].upper() for t in terms]
    if len(set(prefixes)) != len(prefixes):
        raise ValueError(
            f"A condition takes at most one term per ontology type, got {terms}. "
            "Run separate queries to compare tissues."
        )
    return "_".join(terms)


def variant_window(variant: GenomicVariant, length: int) -> GenomicInterval:
    """`length`-bp window centred on the REF allele, as `variant.reference_interval.resize`."""
    return variant.reference_interval.resized(length)


# --------------------------------------------------------------------------- cache


class ResponseCache:
    """Tiny content-addressed JSON cache: one file per request payload."""

    def __init__(self, directory: Path, ttl_seconds: float | None = None) -> None:
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.ttl_seconds = ttl_seconds

    @staticmethod
    def key_for(payload: dict[str, Any]) -> str:
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        return hashlib.sha256(canonical.encode()).hexdigest()

    def _path(self, key: str) -> Path:
        return self.directory / f"{key}.json"

    def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:
        path = self._path(self.key_for(payload))
        if not path.exists():
            return None
        if self.ttl_seconds is not None and time.time() - path.stat().st_mtime > self.ttl_seconds:
            logger.debug("Cache entry expired: %s", path.name)
            return None
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            logger.warning("Corrupt cache entry %s, ignoring", path.name)
            return None

    def set(self, payload: dict[str, Any], value: dict[str, Any]) -> Path:
        path = self._path(self.key_for(payload))
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(value, indent=2, default=str))
        tmp.replace(path)  # atomic on POSIX
        return path


# --------------------------------------------------------------------------- client


def _import_sdk() -> tuple[Any, Any]:
    try:
        from alphagenome.data import genome
        from alphagenome.models.v2 import dna_client
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise AlphaGenomeUnavailableError(
            "The AlphaGenome SDK is not installed. Run: pip install -e '.[alphagenome]'"
        ) from exc
    return genome, dna_client


class AlphaGenomeClient:
    """Cached wrapper around `alphagenome.models.dna_client`."""

    def __init__(
        self,
        api_key: str,
        cache_dir: Path,
        *,
        sequence_length: int = 1_048_576,
        use_cache: bool = True,
        min_interval_s: float = 1.0,
        cache_ttl_seconds: float | None = None,
    ) -> None:
        if not api_key:
            raise ValueError("api_key must be a non-empty string")
        self._validate_length(sequence_length)
        self._api_key = api_key  # never logged
        self.sequence_length = sequence_length
        self.use_cache = use_cache
        self.min_interval_s = min_interval_s
        self.cache = ResponseCache(Path(cache_dir) / "alphagenome", ttl_seconds=cache_ttl_seconds)
        self._model: Any | None = None
        self._last_call = 0.0

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(sequence_length={self.sequence_length}, "
            f"use_cache={self.use_cache}, api_key=***)"
        )

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _validate_length(length: int) -> None:
        if length not in SUPPORTED_SEQUENCE_LENGTHS:
            raise ValueError(
                f"sequence_length {length} not supported; "
                f"choose one of {SUPPORTED_SEQUENCE_LENGTHS}"
            )

    @staticmethod
    def _validate_outputs(output_types: Sequence[str]) -> list[str]:
        names = [o.upper() for o in output_types]
        bad = sorted(set(names) - set(TRACK_OUTPUT_TYPES))
        if bad:
            raise ValueError(f"Unsupported output types {bad}; choose from {TRACK_OUTPUT_TYPES}")
        return names

    @property
    def model(self) -> Any:
        if self._model is None:
            _, dna_client = _import_sdk()
            logger.info("Connecting to AlphaGenome API")
            self._model = dna_client.create(self._api_key)
        return self._model

    def _throttle(self) -> None:
        wait = self.min_interval_s - (time.monotonic() - self._last_call)
        if wait > 0:
            logger.debug("Throttling AlphaGenome call for %.2fs", wait)
            time.sleep(wait)
        self._last_call = time.monotonic()

    def _cached_or_call(
        self,
        payload: dict[str, Any],
        call: Callable[[], ModelT],
        model: type[BaseModel] = AlphaGenomePrediction,
    ) -> ModelT:
        if self.use_cache and (hit := self.cache.get(payload)) is not None:
            logger.info("AlphaGenome cache hit (%s) for %s", payload["kind"], payload["query"])
            cached = model.model_validate(hit).model_copy(update={"from_cache": True})
            return cached  # type: ignore[return-value]

        self._throttle()
        logger.info("Querying AlphaGenome (%s) for %s", payload["kind"], payload["query"])
        prediction = call()
        if self.use_cache:
            path = self.cache.set(payload, prediction.model_dump(mode="json"))
            logger.debug("Cached AlphaGenome response at %s", path)
        return prediction

    # ------------------------------------------------------------------ public API

    def predict_variant(
        self,
        variant: GenomicVariant,
        ontology_terms: Sequence[str],
        output_types: Sequence[str] = ("RNA_SEQ",),
        sequence_length: int | None = None,
    ) -> AlphaGenomePrediction:
        """Predict REF vs ALT tracks in a window centred on `variant`."""
        length = sequence_length or self.sequence_length
        self._validate_length(length)
        outputs = self._validate_outputs(output_types)
        to_condition(ontology_terms)  # fail fast, before any network call
        interval = variant_window(variant, length)
        payload = {
            "v": CACHE_SCHEMA_VERSION,
            "kind": "predict_variant",
            "query": str(variant),
            "interval": str(interval),
            "ontology_terms": sorted(ontology_terms),
            "outputs": sorted(outputs),
        }

        def call() -> AlphaGenomePrediction:
            genome, dna_client = _import_sdk()
            result = self.model.predict_variant(
                interval=genome.Interval(interval.chromosome, interval.start, interval.end),
                variant=genome.Variant(
                    chromosome=variant.chromosome,
                    position=variant.position,
                    reference_bases=variant.reference_bases,
                    alternate_bases=variant.alternate_bases,
                ),
                condition=to_condition(ontology_terms),
                requested_outputs=[getattr(dna_client.OutputType, o) for o in outputs],
            )
            tracks = [
                s
                for o in outputs
                if (
                    s := summarize_variant_tracks(
                        o,
                        getattr(result.reference, o.lower(), None),
                        getattr(result.alternate, o.lower(), None),
                    )
                )
                is not None
            ]
            return AlphaGenomePrediction(
                query=str(variant),
                interval=interval,
                ontology_terms=list(ontology_terms),
                tracks=tracks,
            )

        return self._cached_or_call(payload, call)

    def predict_interval(
        self,
        interval: GenomicInterval,
        ontology_terms: Sequence[str],
        output_types: Sequence[str] = ("RNA_SEQ",),
        sequence_length: int | None = None,
    ) -> AlphaGenomePrediction:
        """Predict tracks over `interval`, resized (centred) to a supported length."""
        length = sequence_length or self.sequence_length
        self._validate_length(length)
        if interval.width > length:
            logger.warning(
                "Interval %s (%d bp) is wider than the %d bp context; it will be truncated",
                interval,
                interval.width,
                length,
            )
        outputs = self._validate_outputs(output_types)
        to_condition(ontology_terms)  # fail fast, before any network call
        window = interval.resized(length)
        payload = {
            "v": CACHE_SCHEMA_VERSION,
            "kind": "predict_interval",
            "query": str(interval),
            "interval": str(window),
            "ontology_terms": sorted(ontology_terms),
            "outputs": sorted(outputs),
        }

        def call() -> AlphaGenomePrediction:
            genome, dna_client = _import_sdk()
            result = self.model.predict_interval(
                interval=genome.Interval(window.chromosome, window.start, window.end),
                condition=to_condition(ontology_terms),
                requested_outputs=[getattr(dna_client.OutputType, o) for o in outputs],
            )
            tracks = [
                s
                for o in outputs
                if (s := summarize_interval_tracks(o, getattr(result, o.lower(), None))) is not None
            ]
            return AlphaGenomePrediction(
                query=str(interval),
                interval=window,
                ontology_terms=list(ontology_terms),
                tracks=tracks,
            )

        return self._cached_or_call(payload, call)

    def score_variant(
        self,
        variant: GenomicVariant,
        ontology_terms: Sequence[str],
        sequence_length: int | None = None,
    ) -> VariantScores:
        """Score `variant` with AlphaGenome's recommended variant scorers.

        Unlike `predict_variant`, this uses the model's own scoring recipes
        (centre-masked log2 ratios, gene-level LFC, splicing scores, ...).
        These are the numbers to report and compare across variants.
        """
        length = sequence_length or self.sequence_length
        self._validate_length(length)
        condition = to_condition(ontology_terms)
        interval = variant_window(variant, length)
        payload = {
            "v": CACHE_SCHEMA_VERSION,
            "kind": "score_variant",
            "query": str(variant),
            "interval": str(interval),
            "condition": condition,
        }

        def call() -> VariantScores:
            genome, _ = _import_sdk()
            from alphagenome.models.v2 import variant_scorers

            scores = self.model.score_variant(
                genome.Interval(interval.chromosome, interval.start, interval.end),
                genome.Variant(
                    chromosome=variant.chromosome,
                    position=variant.position,
                    reference_bases=variant.reference_bases,
                    alternate_bases=variant.alternate_bases,
                ),
                condition=condition,
            )
            records = records_from_tidy(variant_scorers.tidy_scores(scores))
            return VariantScores(
                query=str(variant),
                interval=interval,
                condition=condition,
                summary=summarize_scores(records),
                records=records,
            )

        return self._cached_or_call(payload, call, VariantScores)


# --------------------------------------------------------------------------- reducers


def _track_names(track_data: Any, n: int) -> list[str]:
    meta = getattr(track_data, "metadata", None)
    if meta is not None and "name" in getattr(meta, "columns", []):
        return [str(x) for x in meta["name"].tolist()]
    return [f"track_{i}" for i in range(n)]


def summarize_variant_tracks(output_type: str, ref: Any, alt: Any) -> TrackSummary | None:
    """Collapse REF/ALT `TrackData` (positions x tracks) into a `TrackSummary`."""
    if ref is None or alt is None:
        logger.warning("No %s tracks returned (check ontology terms)", output_type)
        return None
    ref_v = np.asarray(ref.values, dtype=float)
    alt_v = np.asarray(alt.values, dtype=float)
    if ref_v.size == 0:
        return None
    delta = alt_v - ref_v
    per_track_max = np.abs(delta).max(axis=0)
    top = int(per_track_max.argmax())
    ref_sum, alt_sum = float(ref_v.sum()), float(alt_v.sum())
    return TrackSummary(
        output_type=output_type,
        n_tracks=int(ref_v.shape[1]) if ref_v.ndim > 1 else 1,
        resolution_bp=getattr(ref, "resolution", None),
        ref_sum=ref_sum,
        alt_sum=alt_sum,
        log2_fold_change=float(np.log2((alt_sum + 1.0) / (ref_sum + 1.0))),
        max_abs_delta=float(per_track_max[top]),
        top_track=_track_names(ref, per_track_max.size)[top],
    )


def summarize_interval_tracks(output_type: str, data: Any) -> TrackSummary | None:
    if data is None:
        logger.warning("No %s tracks returned (check ontology terms)", output_type)
        return None
    values = np.asarray(data.values, dtype=float)
    if values.size == 0:
        return None
    per_track_mean = values.mean(axis=0)
    top = int(np.atleast_1d(per_track_mean).argmax())
    return TrackSummary(
        output_type=output_type,
        n_tracks=int(values.shape[1]) if values.ndim > 1 else 1,
        resolution_bp=getattr(data, "resolution", None),
        mean_signal=float(values.mean()),
        top_track=_track_names(data, np.atleast_1d(per_track_mean).size)[top],
    )


# --------------------------------------------------------------------------- variant scores

_SCORER_NOISE = re.compile(r"requested_output=\w+,?\s*|aggregation_type=")


def short_scorer_name(scorer: str) -> str:
    """Drop redundant fields from a scorer repr.

    'CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=ACTIVE_SUM)'
    -> 'CenterMaskScorer(width=501, ACTIVE_SUM)'
    """
    name = _SCORER_NOISE.sub("", scorer)
    return name.replace(", )", ")").replace("(, ", "(")


def scorer_kind(scorer: str) -> str:
    """ACTIVE_* / *ActiveScorer measure activity (max of REF, ALT), not ALT-vs-REF effect."""
    return "activity" if "Active" in scorer or "ACTIVE" in scorer else "effect"


def records_from_tidy(tidy: pd.DataFrame | None) -> list[ScoreRecord]:
    """Keep the columns we report from `variant_scorers.tidy_scores` output."""
    if tidy is None or tidy.empty:
        return []
    text_cols = ["track_name", "gene_name", "biosample_name"]
    df = tidy.reindex(columns=["output_type", "variant_scorer", *text_cols, "raw_score"])
    df = df.dropna(subset=["raw_score"]).copy()
    df["output_type"] = df["output_type"].astype(str)
    df["variant_scorer"] = df["variant_scorer"].astype(str).map(short_scorer_name)
    df[text_cols] = df[text_cols].astype(object).where(df[text_cols].notna(), None)
    return [ScoreRecord(**row) for row in df.to_dict(orient="records")]


def _text_or_none(value: Any) -> str | None:
    return None if value is None or pd.isna(value) else str(value)


def summarize_scores(records: Sequence[ScoreRecord]) -> list[ScorerSummary]:
    """Strongest |score| per (output type, scorer); effect scorers listed first."""
    if not records:
        return []
    df = pd.DataFrame([r.model_dump() for r in records])
    df["abs_score"] = df["raw_score"].abs()
    out = []
    for (otype, scorer), group in df.groupby(["output_type", "variant_scorer"], sort=True):
        top = group.loc[group["abs_score"].idxmax()]
        out.append(
            ScorerSummary(
                output_type=str(otype),
                scorer=str(scorer),
                kind=scorer_kind(str(scorer)),
                n_scores=len(group),
                max_abs_score=float(top["abs_score"]),
                top_score=float(top["raw_score"]),
                top_track=_text_or_none(top["track_name"]),
                top_gene=_text_or_none(top["gene_name"]),
            )
        )
    return sorted(out, key=lambda x: (x.kind != "effect", x.output_type, x.scorer))
