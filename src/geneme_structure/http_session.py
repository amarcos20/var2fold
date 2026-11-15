"""Shared `requests` session with retries and a polite User-Agent."""

from __future__ import annotations

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from geneme_structure import __version__

USER_AGENT = f"geneme-structure/{__version__} (python-requests)"


def build_session(total_retries: int = 4, backoff: float = 1.0) -> requests.Session:
    """Session that retries 429/5xx with exponential backoff (honours Retry-After)."""
    retry = Retry(
        total=total_retries,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        respect_retry_after_header=True,
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    return session
