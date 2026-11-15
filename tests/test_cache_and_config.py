from types import SimpleNamespace

import numpy as np
import pytest

from var2fold.alphagenome_client import (
    AlphaGenomeClient,
    ResponseCache,
    summarize_variant_tracks,
)
from var2fold.config import MissingCredentialError, Settings


def test_cache_roundtrip_is_order_independent(tmp_path):
    cache = ResponseCache(tmp_path)
    cache.set({"a": 1, "b": [1, 2]}, {"ok": True})
    assert cache.get({"b": [1, 2], "a": 1}) == {"ok": True}
    assert cache.get({"a": 2}) is None


def test_placeholder_key_is_treated_as_missing(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)  # no .env here
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "your_api_key_here")
    s = Settings.from_env()
    assert not s.has_alphagenome_key
    with pytest.raises(MissingCredentialError):
        s.require_alphagenome_key()


def test_key_never_in_repr(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "super-secret")
    s = Settings.from_env()
    assert "super-secret" not in repr(s)
    client = AlphaGenomeClient(s.require_alphagenome_key(), tmp_path)
    assert "super-secret" not in repr(client)


def test_cache_hit_skips_api(tmp_path):
    from var2fold.models import GenomicVariant

    client = AlphaGenomeClient("k", tmp_path, sequence_length=16_384, min_interval_s=0)
    calls = []

    def fake_call():
        from var2fold.models import AlphaGenomePrediction, GenomicInterval

        calls.append(1)
        return AlphaGenomePrediction(
            query="q",
            interval=GenomicInterval(chromosome="chr1", start=0, end=16_384),
            condition="X",
            tracks=[],
        )

    payload = {"kind": "test", "query": str(GenomicVariant.parse("chr1:100:A>G"))}
    first = client._cached_or_call(payload, fake_call)
    second = client._cached_or_call(payload, fake_call)
    assert len(calls) == 1 and not first.from_cache and second.from_cache


def test_summarize_variant_tracks():
    import pandas as pd

    ref = SimpleNamespace(
        values=np.zeros((4, 2)), resolution=1, metadata=pd.DataFrame({"name": ["t0", "t1"]})
    )
    alt = SimpleNamespace(values=np.array([[0, 0], [0, 3], [0, 0], [1, 0]], float))
    s = summarize_variant_tracks("RNA_SEQ", ref, alt)
    assert s.n_tracks == 2 and s.top_track == "t1" and s.max_abs_delta == 3.0
    assert s.log2_fold_change == pytest.approx(np.log2(5.0))


def test_to_condition():
    from var2fold.alphagenome_client import to_condition

    expected = "UBERON:0000178_CL:0001054"
    assert to_condition(["UBERON:0000178", "CL:0001054"]) == expected
    assert to_condition("UBERON:0000178_CL:0001054") == expected
    assert to_condition("UBERON:0000178 CL:0001054") == expected
    with pytest.raises(ValueError):
        to_condition("colon")
    with pytest.raises(ValueError):
        to_condition(["UBERON:0000178", "UBERON:0001157"])
    with pytest.raises(ValueError):
        to_condition([])
