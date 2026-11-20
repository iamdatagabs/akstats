"""Public data download, parsing, and schema validation."""
from __future__ import annotations
import gzip
import logging
from io import StringIO
from pathlib import Path

import pandas as pd
import requests

from akstats.config import (
    GEO_MATRIX_URL, PANEL_HEADER_ROW, PANEL_SHEET, RAW_DIR, VARIANTS_HEADER_ROW,
    VARIANTS_SHEET, PANEL_URL, VARIANTS_URL,
)

LOGGER = logging.getLogger(__name__)

# Original supplementary workbooks committed to data/raw/.
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
    required_variants = {"Case", "Hugo_Symbol", "Variant_Classification",
                         "Normal_MAF", "AK_MAF", "SCC_MAF"}
    missing_panel = required_panel - set(panel.columns)
    missing_variants = required_variants - set(variants.columns)
    if missing_panel or missing_variants:
        raise ValueError(f"Source schema mismatch: panel={sorted(missing_panel)}, variants={sorted(missing_variants)}")
    for column in ("Normal_MAF", "AK_MAF", "SCC_MAF"):
        variants[column] = pd.to_numeric(variants[column], errors="coerce")
        values = variants[column].dropna()
        if values.lt(0).any() or values.gt(1).any():
            raise ValueError(f"{column} contains values outside [0, 1].")


def load_mutation_sources() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load committed NCBI/article supplementary files, downloading only as fallback."""
    panel_source = RAW_DIR / PANEL_SOURCE_FILENAME
    variants_source = RAW_DIR / VARIANTS_SOURCE_FILENAME

    # Prefer the exact, version-controlled source workbooks. The URL/cache path
    # remains a fallback for fresh clones that do not include the raw assets.
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
    variants = pd.read_excel(variants_path, sheet_name=VARIANTS_SHEET, header=VARIANTS_HEADER_ROW)
    panel = panel.dropna(subset=["Gene"]).copy()
    variants = variants.dropna(subset=["Case", "Hugo_Symbol"]).copy()
    validate_mutation_columns(panel, variants)
    return panel, variants


def load_geo_series_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse GEO metadata and probe-level expression matrix; no gene mapping is assumed."""
    path = download_file(GEO_MATRIX_URL, RAW_DIR / "GSE158634_series_matrix.txt.gz")
    metadata_lines, matrix_lines, in_matrix = [], [], False
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_matrix = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_matrix and not line.startswith("!"):
                matrix_lines.append(line)
            elif not in_matrix and line.startswith((
                "!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1",
                "!Sample_characteristics_ch1",
            )):
                metadata_lines.append(line)
    if not matrix_lines:
        raise ValueError("GEO series matrix has no expression table.")
    parsed: dict[str, list[str]] = {}
    for line in metadata_lines:
        fields = line.rstrip("\n").split("\t")
        key = fields[0].lstrip("!")
        values = [value.strip('"') for value in fields[1:]]
        if key in parsed:
            parsed[key] = [f"{old}; {new}" for old, new in zip(parsed[key], values)]
        else:
            parsed[key] = values
    metadata = pd.DataFrame(parsed)
    if "Sample_geo_accession" not in metadata:
        raise ValueError("GEO metadata is missing Sample_geo_accession.")
    metadata = metadata.set_index("Sample_geo_accession")
    metadata["Clinical_group"] = metadata.apply(_clinical_group, axis=1)
    expression = pd.read_csv(StringIO("".join(matrix_lines)), sep="\t")
    if "ID_REF" not in expression:
        raise ValueError("GEO expression matrix is missing ID_REF.")
    expression = expression.set_index("ID_REF").apply(pd.to_numeric, errors="coerce")
    common = metadata.index.intersection(expression.columns)
    if common.empty:
        raise ValueError("No sample IDs overlap between GEO metadata and expression matrix.")
    return metadata.loc[common], expression.loc[:, common]


def _clinical_group(row: pd.Series) -> str:
    text = " ".join(str(value).lower() for value in row.astype(str).values)
    if "normal" in text:
        return "Normal skin"
    if "actinic" in text:
        return "Actinic keratosis"
    if "squamous" in text or "scc" in text:
        return "Cutaneous SCC"
    return "Unclassified"
