"""Runtime configuration, read strictly from environment variables.

`.env` is loaded with python-dotenv (searched from the current working
directory upwards) but never overrides variables already set in the
environment. API keys are wrapped in `SecretStr` so they never show up in
logs, reprs or tracebacks.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from dotenv import find_dotenv, load_dotenv
from pydantic import BaseModel, ConfigDict, Field, SecretStr

logger = logging.getLogger(__name__)

_PLACEHOLDER_KEYS = {"", "your_api_key_here"}


class MissingCredentialError(RuntimeError):
    """Raised when a required credential is not configured."""


class Settings(BaseModel):
    model_config = ConfigDict(frozen=True)

    alphagenome_api_key: SecretStr | None = None
    ncbi_api_key: SecretStr | None = None  # optional: raises the NCBI E-utilities rate limit
    cache_dir: Path = Path(".cache")
    data_dir: Path = Path("data")
    http_timeout: float = Field(default=30.0, gt=0)

    @classmethod
    def from_env(cls) -> Settings:
        dotenv_path = find_dotenv(usecwd=True)
        if dotenv_path:
            load_dotenv(dotenv_path, override=False)
            logger.debug("Loaded environment from %s", dotenv_path)

        raw_key = os.getenv("ALPHAGENOME_API_KEY")
        key = None if raw_key is None or raw_key.strip() in _PLACEHOLDER_KEYS else raw_key.strip()

        ncbi_key = (os.getenv("NCBI_API_KEY") or "").strip()

        return cls(
            alphagenome_api_key=SecretStr(key) if key else None,
            ncbi_api_key=SecretStr(ncbi_key) if ncbi_key else None,
            cache_dir=Path(os.getenv("GENEME_CACHE_DIR", ".cache")),
            data_dir=Path(os.getenv("GENEME_DATA_DIR", "data")),
            http_timeout=float(os.getenv("GENEME_HTTP_TIMEOUT", "30")),
        )

    @property
    def has_alphagenome_key(self) -> bool:
        return self.alphagenome_api_key is not None

    def require_alphagenome_key(self) -> str:
        if self.alphagenome_api_key is None:
            raise MissingCredentialError(
                "ALPHAGENOME_API_KEY is not set. Copy .env.example to .env and add your key "
                "(get one at https://deepmind.google.com/science/alphagenome)."
            )
        return self.alphagenome_api_key.get_secret_value()
