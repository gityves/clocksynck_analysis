"""Rule-version stamp of the standard tables (common.data): round trip, stale refusal, legacy fallback warning."""
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

sys.path[:0] = [str(Path(__file__).resolve().parents[2]), str(Path(__file__).resolve().parents[1] / "datasets" / "_atlas")]
from common import data  # noqa: E402

NAME = "kidney_sc"


@pytest.fixture
def store(tmp_path, monkeypatch):
    monkeypatch.setattr(data, "PSEUDOBULK", tmp_path)
    pb = pd.DataFrame({"donor": ["a", "b"], "compartment": ["x", "x"], "n_cells": [1, 2]})
    donors = pd.DataFrame({"age": [1, 2]}, index=pd.Index(["a", "b"], name="donor"))
    return pb, donors


def test_roundtrip_is_stamped(store):
    data.save(NAME, *store)
    assert data._stored_hash(data.PSEUDOBULK / f"{NAME}.parquet") == data.rules_hash(NAME)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        pb, donors = data.load(NAME)
    assert list(donors.index) == ["a", "b"]


def test_changed_rules_refuse(store, monkeypatch):
    data.save(NAME, *store)
    monkeypatch.setattr(data, "rules_hash", lambda name: "r1-changed")
    with pytest.raises(data.StaleTableError):
        data.load(NAME)
    with pytest.warns(UserWarning, match="rebuild"):
        data.load(NAME, allow_stale=True)


def test_unstamped_table_warns(store):
    store[0].to_parquet(data.PSEUDOBULK / f"{NAME}.parquet", index=False)
    with pytest.warns(UserWarning, match="no rule stamp"):
        data.load(NAME)


def test_legacy_fallback_warns(store, monkeypatch):
    class Mod:
        @staticmethod
        def from_legacy(**kw):
            return store

    monkeypatch.setattr(data, "_builder", lambda name: Mod)
    with pytest.warns(UserWarning, match="legacy"):
        data.load(NAME)
