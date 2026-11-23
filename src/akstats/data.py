"""Public mutation-data download, parsing, and schema validation."""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import requests

from akstats.config import (
    PANEL_HEADER_ROW,
    PANEL_SHEET,
    PANEL_URL,
    RAW_DIR,
    VARIANTS_HEADER_ROW,
    VARIANTS_SHEET,
    VARIANTS_URL,
)

LOGGER = logging.getLogger(__name__)

PANEL_SOURCE_FILENAME = "41467_2025_65687_MOESM8_ESM.xlsx"
VARIANTS_SOURCE_FILENAME = "41467_2025_65687_MOESM9_ESM.xlsx"


def download_file(url: str, destination: Path, timeout: int = 90) -> Path:
    """Download atomically and reuse a non-empty local cache."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size:
        return destination
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    temporary = destination.with_suffix(destination.suffix + ".part")
    temporary.write_bytes(response.content)
    temporary.replace(destination)
    LOGGER.info("Downloaded %s", destination.name)
    return destination


def validate_mutation_columns(panel: pd.DataFrame, variants: pd.DataFrame) -> None:
    required_panel = {"Gene"}
    required_variants = {
        "Case", "Hugo_Symbol", "Variant_Classification",
        "Normal_MAF", "AK_MAF", "SCC_MAF",
    }
    missing_panel = required_panel - set(panel.columns)
    missing_variants = required_variants - set(variants.columns)
    if missing_panel or missing_variants:
        raise ValueError(
            "Source schema mismatch: "
            f"panel={sorted(missing_panel)}, variants={sorted(missing_variants)}"
        )
    for column in ("Normal_MAF", "AK_MAF", "SCC_MAF"):
        variants[column] = pd.to_numeric(variants[column], errors="coerce")
        values = variants[column].dropna()
        if values.lt(0).any() or values.gt(1).any():
            raise ValueError(f"{column} contains values outside [0, 1].")


def load_mutation_sources() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load committed supplementary files, downloading only as fallback."""
    panel_source = RAW_DIR / PANEL_SOURCE_FILENAME
    variants_source = RAW_DIR / VARIANTS_SOURCE_FILENAME
    panel_path = (
        panel_source if panel_source.is_file()
        else download_file(PANEL_URL, RAW_DIR / "gene_panel.xlsx")
    )
    variants_path = (
        variants_source if variants_source.is_file()
        else download_file(VARIANTS_URL, RAW_DIR / "somatic_variants.xlsx")
    )
    LOGGER.info("Panel source: %s", panel_path.name)
    LOGGER.info("Variant source: %s", variants_path.name)
    panel = pd.read_excel(panel_path, sheet_name=PANEL_SHEET, header=PANEL_HEADER_ROW)
    variants = pd.read_excel(
        variants_path, sheet_name=VARIANTS_SHEET, header=VARIANTS_HEADER_ROW
    )
    panel = panel.dropna(subset=["Gene"]).copy()
    variants = variants.dropna(subset=["Case", "Hugo_Symbol"]).copy()
    validate_mutation_columns(panel, variants)
    return panel, variants
