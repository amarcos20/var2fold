"""Shared pytest fixtures."""

import pytest
from fakes import FakeSession


@pytest.fixture
def fake_session() -> FakeSession:
    return FakeSession()
