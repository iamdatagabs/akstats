"""Independent exploratory transcriptomics layer for GEO GSE158634."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from akstats.config import GENE_TO_MODULE, RAW_DIR

LOGGER = logging.getLogger(__name__)
GEO_COUNT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158634/"
    "suppl/GSE158634_3.count.csv.gz"
)
GEO_COUNT_FILE = "GSE158634_3.count.csv.gz"
SAMPLES = {
    "GSM4804794": ("NS_001", "Normal skin"),
    "GSM4804795": ("NS_002", "Normal skin"),
    "GSM4804796": ("NS_004", "Normal skin"),
    "GSM4804797": ("NS_005", "Normal skin"),
    "GSM4804798": ("NS_006", "Normal skin"),
    "GSM4804799": ("NS_007", "Normal skin"),
    "GSM4804800": ("NS_008", "Normal skin"),
    "GSM4804801": ("NS_009", "Normal skin"),
    "GSM4804802": ("AK_010", "Actinic keratosis"),
    "GSM4804803": ("AK_011", "Actinic keratosis"),
    "GSM4804804": ("AK_012", "Actinic keratosis"),
    "GSM4804805": ("AK_013", "Actinic keratosis"),
    "GSM4804806": ("AK_014", "Actinic keratosis"),
    "GSM4804807": ("AK_015", "Actinic keratosis"),
    "GSM4804808": ("AK_016", "Actinic keratosis"),
    "GSM4804809": ("AK_017", "Actinic keratosis"),
    "GSM4804810": ("cSCC_018", "Cutaneous SCC"),
    "GSM4804811": ("cSCC_019", "Cutaneous SCC"),
    "GSM4804812": ("cSCC_020", "Cutaneous SCC"),
    "GSM4804813": ("cSCC_021", "Cutaneous SCC"),
}
TARGETS = {gene[0] + gene[1:].lower(): module for gene, module in GENE_TO_MODULE.items()}
TARGETS["Rora"] = "Literature candidate"


def _download_counts(destination: Path) -> Path:
    import requests

    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists() or not destination.stat().st_size:
        response = requests.get(GEO_COUNT_URL, timeout=120)
        response.raise_for_status()
        temporary = destination.with_suffix(destination.suffix + ".part")
        temporary.write_bytes(response.content)
        temporary.replace(destination)
    return destination


def load_geo_counts(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load integer gene counts and map the 20 GEO sample accessions to tissue stages."""
    matrix = pd.read_csv(path, compression="gzip")
    if matrix.shape[1] < 21:
        raise ValueError(f"GEO count matrix has only {matrix.shape[1] - 1} sample columns.")
    feature_col = matrix.columns[0]
    matrix = matrix.rename(columns={feature_col: "feature_id"}).dropna(subset=["feature_id"])
    matrix["feature_id"] = matrix["feature_id"].astype(str).str.strip()
    matrix = matrix.set_index("feature_id")
    selected = {}
    metadata = []
    for accession, (title, stage) in SAMPLES.items():
        candidates = [c for c in matrix.columns if str(c).strip() in {accession, title}]
        if not candidates:
            candidates = [c for c in matrix.columns if title in str(c) or accession in str(c)]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one count column for {accession}/{title}; found {len(candidates)}."
            )
        selected[accession] = candidates[0]
        metadata.append({
            "sample_accession": accession,
            "sample_title": title,
            "stage": stage,
            "animal_id": None,
            "experimental_unit_note": "Animal ID not recoverable from deposited sample labels",
        })
    counts = matrix.loc[:, list(selected.values())].copy()
    counts.columns = list(selected.keys())
    counts = counts.apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any() or (counts < 0).any().any():
        raise ValueError("GEO count matrix contains missing, non-numeric or negative values.")
    if not np.allclose(counts.to_numpy(), np.round(counts.to_numpy())):
        raise ValueError("Supplementary matrix is not integer raw counts; count model not applied.")
    counts = counts.round().astype("int64")
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    return counts, pd.DataFrame(metadata)


def summarize_geo_counts(counts: pd.DataFrame, metadata: pd.DataFrame) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    """Return sample QC, target-gene log2-CPM summaries and descriptive stage contrasts."""
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("At least one GEO library has zero total counts.")
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 0.5)
    qc = metadata.copy()
    qc["library_size"] = qc["sample_accession"].map(library_sizes).astype(int)
    qc["detected_features"] = qc["sample_accession"].map((counts > 0).sum()).astype(int)

    canonical = {str(gene).casefold(): gene for gene in TARGETS}
    matches = {feature: canonical[feature.casefold()] for feature in log_cpm.index
               if feature.casefold() in canonical}
    if not matches:
        return (
            qc,
            pd.DataFrame(columns=["Gene", "Module", "Stage", "n_samples", "median_log2_cpm", "q1_log2_cpm", "q3_log2_cpm"]),
            pd.DataFrame(columns=["Gene", "Module", "contrast", "median_log2_cpm_difference", "interpretation"]),
        )
    rows = []
    for feature, gene in matches.items():
        module = TARGETS[gene]
        for stage in ("Normal skin", "Actinic keratosis", "Cutaneous SCC"):
            samples = metadata.loc[metadata["stage"].eq(stage), "sample_accession"].tolist()
            values = log_cpm.loc[feature, samples].astype(float)
            rows.append({
                "Gene": gene,
                "Module": module,
                "Stage": stage,
                "n_samples": len(samples),
                "median_log2_cpm": float(values.median()),
                "q1_log2_cpm": float(values.quantile(.25)),
                "q3_log2_cpm": float(values.quantile(.75)),
            })
    summary = pd.DataFrame(rows)
    contrasts = []
    for gene, frame in summary.groupby("Gene"):
        values = frame.set_index("Stage")["median_log2_cpm"]
        for left, right, label in (
            ("Normal skin", "Actinic keratosis", "AK_vs_Normal"),
            ("Actinic keratosis", "Cutaneous SCC", "cSCC_vs_AK"),
            ("Normal skin", "Cutaneous SCC", "cSCC_vs_Normal"),
        ):
            if left in values and right in values:
                contrasts.append({
                    "Gene": gene,
                    "Module": frame["Module"].iloc[0],
                    "contrast": label,
                    "median_log2_cpm_difference": float(values[right] - values[left]),
                    "interpretation": "Descriptive difference of group medians; not differential-expression inference",
                })
    return qc, summary, pd.DataFrame(contrasts)


def run_transcriptomics(output_dir: Path) -> dict:
    """Run the independent GEO QC and descriptive target-expression workflow."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "geo_sample_qc.csv").write_text("", encoding="utf-8")
    pd.DataFrame(columns=["Gene", "Module", "Stage", "n_samples", "median_log2_cpm", "q1_log2_cpm", "q3_log2_cpm"]).to_csv(output_dir / "geo_target_expression.csv", index=False)
    pd.DataFrame(columns=["Gene", "Module", "contrast", "median_log2_cpm_difference", "interpretation"]).to_csv(output_dir / "geo_target_contrasts.csv", index=False)
    status = {
        "accession": "GSE158634",
        "organism": "Mus musculus",
        "assay": "bulk RNA-seq",
        "status": "unavailable",
        "n_samples_expected": len(SAMPLES),
        "animal_id_available": False,
        "inference_status": "not_performed",
        "notice": (
            "Independent mouse cohort. Sample-level summaries are descriptive only; "
            "animal identifiers are unavailable, so no inferential differential-expression "
            "or enrichment tests are reported."
        ),
    }
    try:
        path = _download_counts(RAW_DIR / GEO_COUNT_FILE)
        counts, metadata = load_geo_counts(path)
        qc, summary, contrasts = summarize_geo_counts(counts, metadata)
        qc.to_csv(output_dir / "geo_sample_qc.csv", index=False)
        summary.to_csv(output_dir / "geo_target_expression.csv", index=False)
        contrasts.to_csv(output_dir / "geo_target_contrasts.csv", index=False)
        status.update({
            "status": "validated" if len(summary) else "counts_validated_gene_ids_unmapped",
            "n_samples": int(counts.shape[1]),
            "n_features": int(counts.shape[0]),
            "identifier_type": "gene_symbol" if len(summary) else "unresolved",
            "n_target_genes": int(summary["Gene"].nunique()) if len(summary) else 0,
            "n_library_qc_rows": len(qc),
        })
    except (OSError, ValueError, KeyError, ImportError) as exc:
        LOGGER.warning("GSE158634 layer was not integrated: %s", exc)
        status["status"] = "validation_failed"
        status["validation_message"] = str(exc)
        pd.DataFrame(columns=["sample_accession", "sample_title", "stage", "animal_id", "experimental_unit_note", "library_size", "detected_features"]).to_csv(output_dir / "geo_sample_qc.csv", index=False)
    (output_dir / "geo_transcriptomics_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return status
