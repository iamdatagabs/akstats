"""NCBI Entrez E-utilities and GEO programmatic-access client."""
from __future__ import annotations

import gzip
import logging
import os
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

LOGGER = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_SERIES_SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158634/soft/"
    "GSE158634_family.soft.gz"
)
GEO_COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158634/suppl/"
    "GSE158634_3.count.csv.gz"
)


class NCBIError(RuntimeError):
    """Raised when NCBI data cannot be retrieved or validated."""


class NCBIGEOClient:
    """Small, cache-aware client for NCBI GEO metadata and supplementary files.

    E-utilities are used for Entrez discovery/summary. GEO's official FTP endpoint
    is used for complete SOFT metadata and large supplementary count matrices,
    as recommended in the GEO programmatic-access documentation.
    """

    def __init__(
        self,
        email: str | None = None,
        api_key: str | None = None,
        tool: str = "akstats",
        timeout: int = 90,
        session: requests.Session | None = None,
    ) -> None:
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.tool = tool
        self.timeout = timeout
        self.session = session or requests.Session()
        retry = Retry(
            total=4,
            connect=4,
            read=4,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))
        self.session.headers.update({"User-Agent": f"{tool}/0.1 (NCBI GEO research client)"})

    def _eutils_params(self, params: dict[str, Any]) -> dict[str, Any]:
        result = {"tool": self.tool, **params}
        if self.email:
            result["email"] = self.email
        if self.api_key:
            result["api_key"] = self.api_key
        return result

    def _get(self, url: str, params: dict[str, Any] | None = None) -> requests.Response:
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            raise NCBIError(f"NCBI request failed for {url}: {exc}") from exc

    def search_geo_series(self, accession: str = "GSE158634") -> list[str]:
        """Search the GEO DataSets Entrez database for a GEO series accession."""
        response = self._get(
            f"{EUTILS_BASE}/esearch.fcgi",
            self._eutils_params({
                "db": "gds",
                "term": f"{accession}[ACCN]",
                "retmode": "json",
                "retmax": 20,
            }),
        )
        try:
            return response.json()["esearchresult"]["idlist"]
        except (ValueError, KeyError, TypeError) as exc:
            raise NCBIError("NCBI ESearch returned an unexpected GEO response.") from exc

    def summarize_geo_series(self, uids: list[str]) -> dict[str, Any]:
        """Retrieve Entrez GEO DataSet summaries for the supplied UIDs."""
        if not uids:
            return {}
        response = self._get(
            f"{EUTILS_BASE}/esummary.fcgi",
            self._eutils_params({
                "db": "gds",
                "id": ",".join(uids),
                "retmode": "json",
            }),
        )
        try:
            return response.json()["result"]
        except (ValueError, KeyError, TypeError) as exc:
            raise NCBIError("NCBI ESummary returned an unexpected GEO response.") from exc

    def download(self, url: str, destination: Path) -> Path:
        """Download a GEO file atomically and reuse a non-empty local cache."""
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.is_file() and destination.stat().st_size:
            return destination
        response = self._get(url)
        temporary = destination.with_suffix(destination.suffix + ".part")
        temporary.write_bytes(response.content)
        if not temporary.stat().st_size:
            temporary.unlink(missing_ok=True)
            raise NCBIError(f"NCBI returned an empty file for {url}")
        temporary.replace(destination)
        LOGGER.info("Downloaded NCBI resource: %s", destination.name)
        return destination

    def fetch_series_soft(self, destination: Path) -> Path:
        """Download the complete GEO family SOFT file (sample characteristics)."""
        return self.download(GEO_SERIES_SOFT_URL, destination)

    def fetch_count_matrix(self, destination: Path) -> Path:
        """Download the deposited integer gene-count matrix."""
        return self.download(GEO_COUNTS_URL, destination)

    def fetch_study(self, cache_dir: Path, accession: str = "GSE158634") -> dict[str, Any]:
        """Fetch Entrez summary, family SOFT metadata and supplementary counts."""
        if accession != "GSE158634":
            raise ValueError("This configured AKStats workflow currently supports GSE158634 only.")
        uids = self.search_geo_series(accession)
        summaries = self.summarize_geo_series(uids)
        soft_path = self.fetch_series_soft(cache_dir / f"{accession}_family.soft.gz")
        counts_path = self.fetch_count_matrix(cache_dir / "GSE158634_3.count.csv.gz")
        return {
            "accession": accession,
            "entrez_uids": uids,
            "entrez_summaries": summaries,
            "soft_path": soft_path,
            "counts_path": counts_path,
        }


def read_gzip_text(path: Path) -> str:
    """Read a UTF-8 gzip text resource with a clear validation error."""
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    except (OSError, EOFError) as exc:
        raise NCBIError(f"Unable to read compressed GEO metadata: {path}") from exc
