"""Design-aware, reproducible transcriptomics workflow for GEO GSE158634."""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from akstats.config import GENE_TO_MODULE, RAW_DIR
from akstats.ncbi import NCBIError, NCBIGEOClient, read_gzip_text

LOGGER = logging.getLogger(__name__)
GEO_COUNT_FILE = "GSE158634_3.count.csv.gz"
GEO_SOFT_FILE = "GSE158634_family.soft.gz"

# Expected accessions are a completeness check, not a source of biological metadata.
EXPECTED_SAMPLES = {
    "GSM4804794", "GSM4804795", "GSM4804796", "GSM4804797",
    "GSM4804798", "GSM4804799", "GSM4804800", "GSM4804801",
    "GSM4804802", "GSM4804803", "GSM4804804", "GSM4804805",
    "GSM4804806", "GSM4804807", "GSM4804808", "GSM4804809",
    "GSM4804810", "GSM4804811", "GSM4804812", "GSM4804813",
}
TARGETS = {gene[0] + gene[1:].lower(): module for gene, module in GENE_TO_MODULE.items()}
TARGETS["Rora"] = "Literature candidate"
SUMMARY_COLUMNS = [
    "Gene", "Module", "Analysis_group", "Stage", "n_samples",
    "median_log2_cpm", "q1_log2_cpm", "q3_log2_cpm",
]
CONTRAST_COLUMNS = [
    "Gene", "Module", "contrast", "median_log2_cpm_difference", "interpretation",
]


def _parse_characteristic(value: str) -> tuple[str, str]:
    if ":" in value:
        key, item = value.split(":", 1)
        return key.strip().casefold(), item.strip()
    return "characteristic", value.strip()


def _stage_from_text(text: str) -> str:
    value = text.casefold()
    if re.search(r"cscc|cutaneous squamous|squamous cell carcinoma|carcinoma", value):
        return "Cutaneous SCC"
    if re.search(r"actinic keratos|\bak[_ -]?", value):
        return "Actinic keratosis"
    if re.search(r"normal|\bns[_ -]?", value):
        return "Normal skin"
    return "Unclassified"


def _exposure_from_characteristics(characteristics: dict[str, list[str]]) -> str | None:
    candidates = set()
    for key, values in characteristics.items():
        if any(token in key for token in ("uv", "treatment", "exposure", "irradiation")):
            value = " ".join(values).casefold()
            if re.search(r"\buvf\b|unirradiated|non.?irradiated|no uv|without uv", value):
                candidates.add("UVF")
            if re.search(r"\buvr\b|irradiated|uv.?exposed|ultraviolet", value):
                candidates.add("UVR")
    return next(iter(candidates)) if len(candidates) == 1 else None


def _animal_from_characteristics(characteristics: dict[str, list[str]]) -> str | None:
    keys = ("animal", "mouse", "subject", "biological replicate", "organism id")
    for key, values in characteristics.items():
        if any(token in key for token in keys):
            for value in values:
                match = re.search(r"(?:mouse|animal|subject)?\s*([A-Za-z]*\d+[A-Za-z0-9_-]*)", value, re.IGNORECASE)
                if match:
                    return match.group(1)
    return None


def parse_geo_family_soft(text: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Parse sample-level metadata from the official GEO family SOFT record."""
    records: list[dict[str, Any]] = []
    series: dict[str, Any] = {}
    current: dict[str, Any] | None = None

    def finish_sample() -> None:
        if current and current.get("sample_accession"):
            characteristics: dict[str, list[str]] = {}
            for raw in current.get("characteristics_raw", []):
                key, value = _parse_characteristic(raw)
                characteristics.setdefault(key, []).append(value)
            title = current.get("sample_title", "")
            searchable = " ".join([title, current.get("source_name", ""), *[
                f"{key}: {'; '.join(values)}" for key, values in characteristics.items()
            ]])
            stage = _stage_from_text(searchable)
            exposure = _exposure_from_characteristics(characteristics)
            animal_id = _animal_from_characteristics(characteristics)
            sra_match = next((
                re.search(r"SRX\d+", relation, re.IGNORECASE)
                for relation in current.get("relations", [])
                if re.search(r"SRX\d+", relation, re.IGNORECASE)
            ), None)
            biosample_match = next((
                re.search(r"BioSample:\s*(SAMN\d+)", relation, re.IGNORECASE)
                for relation in current.get("relations", [])
                if re.search(r"BioSample:\s*(SAMN\d+)", relation, re.IGNORECASE)
            ), None)
            if stage == "Normal skin" and exposure:
                analysis_group = f"Normal skin ({exposure})"
            elif stage == "Normal skin":
                analysis_group = "Normal skin (exposure unknown)"
            else:
                analysis_group = stage
            title_match = re.search(r"\b(NS|AK|cSCC)[_-]?(\d+)\b", title, re.IGNORECASE)
            records.append({
                "sample_accession": current["sample_accession"],
                "biosample_accession": biosample_match.group(1) if biosample_match else None,
                "sra_experiment_accession": sra_match.group(0) if sra_match else None,
                "sra_instrument_model": None,
                "sample_title": title,
                "stage": stage,
                "analysis_group": analysis_group,
                "exposure_group": exposure,
                "animal_id": animal_id,
                "sample_label": title_match.group(0) if title_match else title,
                "source_name": current.get("source_name"),
                "platform_id": current.get("platform_id"),
                "characteristics_json": json.dumps(characteristics, ensure_ascii=False, sort_keys=True),
                "data_processing": " | ".join(current.get("data_processing", [])),
                "treatment_protocol": " | ".join(current.get("treatment_protocol", [])),
                "metadata_source": "NCBI GEO family SOFT",
            })

    for line in text.splitlines():
        if line.startswith("^SERIES ="):
            series["accession"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_title ="):
            series["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_summary ="):
            series.setdefault("summary", []).append(line.split("=", 1)[1].strip())
        elif line.startswith("!Series_overall_design ="):
            series.setdefault("overall_design", []).append(line.split("=", 1)[1].strip())
        elif line.startswith("^SAMPLE ="):
            finish_sample()
            current = {
                "sample_accession": line.split("=", 1)[1].strip(),
                "characteristics_raw": [],
                "data_processing": [],
                "treatment_protocol": [],
                "relations": [],
            }
        elif current is not None and line.startswith("!Sample_"):
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip().strip('"')
            if key == "!Sample_title":
                current["sample_title"] = value
            elif key == "!Sample_source_name_ch1":
                current["source_name"] = value
            elif key == "!Sample_platform_id":
                current["platform_id"] = value
            elif key.startswith("!Sample_characteristics_ch"):
                current["characteristics_raw"].append(value)
            elif key == "!Sample_data_processing":
                current["data_processing"].append(value)
            elif key.startswith("!Sample_treatment_protocol_ch"):
                current["treatment_protocol"].append(value)
            elif key == "!Sample_relation":
                current["relations"].append(value)
    finish_sample()

    metadata = pd.DataFrame(records)
    if metadata.empty:
        raise ValueError("GEO family SOFT contains no parseable sample records.")
    if metadata["sample_accession"].duplicated().any():
        raise ValueError("GEO family SOFT contains duplicate sample accessions.")
    return metadata, series


def merge_sra_metadata(
    metadata: pd.DataFrame, records: dict[str, dict[str, Any]]
) -> pd.DataFrame:
    """Join SRA experiment details to GEO samples using explicit SRX relations."""
    result = metadata.copy()
    for index, row in result.iterrows():
        accession = row.get("sra_experiment_accession")
        record = records.get(accession) if pd.notna(accession) else None
        if not record:
            continue
        if record.get("biosample_accession") and pd.isna(row.get("biosample_accession")):
            result.at[index, "biosample_accession"] = record["biosample_accession"]
        result.at[index, "sra_instrument_model"] = record.get("instrument_model")
        result.at[index, "metadata_source"] = "NCBI GEO SOFT + SRA"
    return result


def merge_biosample_metadata(
    metadata: pd.DataFrame, records: dict[str, dict[str, Any]]
) -> pd.DataFrame:
    """Enrich GEO SOFT rows with explicit attributes from linked NCBI BioSamples."""
    result = metadata.copy()
    result["biosample_attributes_json"] = "{}"
    for index, row in result.iterrows():
        accession = row.get("biosample_accession")
        record = records.get(accession) if pd.notna(accession) else None
        if not record:
            continue
        attributes = record.get("attributes", {})
        result.at[index, "biosample_attributes_json"] = json.dumps(
            attributes, ensure_ascii=False, sort_keys=True
        )
        merged = {}
        try:
            merged.update(json.loads(row.get("characteristics_json") or "{}"))
        except (TypeError, ValueError):
            pass
        for key, values in attributes.items():
            merged.setdefault(key, []).extend(values)
        exposure = _exposure_from_characteristics(merged)
        animal_id = _animal_from_characteristics(merged)
        if exposure:
            result.at[index, "exposure_group"] = exposure
        if animal_id:
            result.at[index, "animal_id"] = animal_id
        if result.at[index, "stage"] == "Normal skin":
            result.at[index, "analysis_group"] = (
                f"Normal skin ({exposure})" if exposure
                else "Normal skin (exposure unknown)"
            )
        result.at[index, "characteristics_json"] = json.dumps(
            merged, ensure_ascii=False, sort_keys=True
        )
        result.at[index, "metadata_source"] = "NCBI GEO family SOFT + BioSample"
    return result


def _metadata_consistency(metadata: pd.DataFrame) -> dict[str, Any]:
    """Surface conflicting genome-build and sequencing-instrument statements."""
    processing = " ".join(metadata.get("data_processing", pd.Series(dtype=str)).dropna().astype(str))
    builds = sorted(set(re.findall(
        r"\b(?:mm\d+|GRCm\d+|hg\d+|GRCh\d+)\b",
        processing,
        re.IGNORECASE,
    )))
    instruments = set(re.findall(
        r"\b(?:HiSeq|NextSeq|NovaSeq)\s*\d+",
        processing,
        re.IGNORECASE,
    ))
    if "sra_instrument_model" in metadata:
        instruments.update(metadata["sra_instrument_model"].dropna().astype(str))
    instruments = sorted(instruments)
    platform_ids = sorted(set(metadata["platform_id"].dropna().astype(str)))
    if len(builds) > 1:
        build_status = "conflicting_source_metadata"
    elif len(builds) == 1:
        build_status = "single_build_mentioned_not_independently_verified"
    else:
        build_status = "not_reported"
    platform_warning = (
        "GEO platform GPL18480 is titled Illumina HiSeq 1500, while sample processing "
        "text mentions HiSeq 2500; verify against the source publication."
        if "GPL18480" in platform_ids and any("2500" in item for item in instruments)
        else None
    )
    return {
        "genome_builds_detected": builds,
        "genome_build": builds[0] if len(builds) == 1 else (
            "conflicting: " + ", ".join(builds) if builds else "not reported"
        ),
        "genome_build_consistency": build_status,
        "instrument_models_reported": instruments,
        "platform_ids": platform_ids,
        "platform_instrument_warning": platform_warning,
    }


def _download_geo_sources(cache_dir: Path, client: NCBIGEOClient | None = None) -> dict[str, Any]:
    client = client or NCBIGEOClient()
    return client.fetch_study(cache_dir)


def load_geo_counts(
    path: Path, metadata: pd.DataFrame | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load integer counts and join columns to metadata by exact GEO accession/title."""
    matrix = pd.read_csv(path, compression="gzip")
    if matrix.shape[1] < 2:
        raise ValueError("GEO count matrix has no sample columns.")
    feature_col = matrix.columns[0]
    matrix = matrix.rename(columns={feature_col: "feature_id"}).dropna(subset=["feature_id"])
    matrix["feature_id"] = matrix["feature_id"].astype(str).str.strip()
    matrix = matrix.set_index("feature_id")
    if metadata is None:
        raise ValueError("GEO SOFT sample metadata is required; sample groups are not hard-coded.")
    selected: dict[str, Any] = {}
    for row in metadata.itertuples(index=False):
        accession = row.sample_accession
        title = row.sample_title
        candidates = [column for column in matrix.columns
                      if str(column).strip() in {accession, title}]
        if not candidates:
            candidates = [column for column in matrix.columns
                          if accession in str(column) or title == str(column).strip()]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one count column for {accession}/{title}; found {len(candidates)}."
            )
        selected[accession] = candidates[0]
    if set(selected) != set(metadata["sample_accession"]):
        raise ValueError("Count matrix and SOFT metadata sample accessions do not match.")
    counts = matrix.loc[:, list(selected.values())].copy()
    counts.columns = list(selected.keys())
    counts = counts.apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any() or (counts < 0).any().any():
        raise ValueError("GEO count matrix contains missing, non-numeric or negative values.")
    if not np.allclose(counts.to_numpy(), np.round(counts.to_numpy())):
        raise ValueError("Supplementary matrix is not integer raw counts.")
    counts = counts.round().astype("int64")
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    return counts, metadata.copy()


def build_publication_candidate_crosswalk(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a provisional GSM-to-publication sample crosswalk for human review.

    The publication reports four UVF controls and four UVR animals, each with
    normal skin, AK-1, AK-2 and cSCC. GEO title numbering follows that same
    block order, but the paper does not publish an explicit GSM-to-subject key.
    Therefore these assignments are candidates only and are never copied into
    the verified animal_id/exposure_group columns used by inferential models.
    """
    result = metadata.copy()
    candidate_columns = {
        "publication_candidate_animal_id": None,
        "publication_candidate_exposure_group": None,
        "publication_candidate_subject": None,
        "publication_candidate_sample_name": None,
        "publication_crosswalk_status": "unmapped",
        "publication_crosswalk_basis": (
            "Provisional alignment of GEO sample-title numbering to Table 2 "
            "sample order in PMID 33664254; not an author-confirmed GSM crosswalk."
        ),
        "publication_crosswalk_review_required": True,
    }
    for column, default in candidate_columns.items():
        result[column] = default

    normal_map = {
        "NS_001": ("UVF_1", "UVF", 1, "NC-1"),
        "NS_002": ("UVF_2", "UVF", 2, "NC-2"),
        "NS_004": ("UVF_3", "UVF", 3, "NC-3"),
        "NS_005": ("UVF_4", "UVF", 4, "NC-4"),
        "NS_006": ("UVR_1", "UVR", 5, "1-NC"),
        "NS_007": ("UVR_2", "UVR", 6, "2-NC"),
        "NS_008": ("UVR_3", "UVR", 7, "3-NC"),
        "NS_009": ("UVR_4", "UVR", 8, "4-NC"),
    }
    for index, row in result.iterrows():
        title = str(row.get("sample_title", "")).strip()
        assignment = normal_map.get(title)
        if assignment:
            animal, exposure, subject, paper_name = assignment
            role = "Normal skin"
        else:
            ak_match = re.fullmatch(r"AK_(\d+)", title, re.IGNORECASE)
            cscc_match = re.fullmatch(r"cSCC_(\d+)", title, re.IGNORECASE)
            if ak_match:
                number = int(ak_match.group(1))
                if 10 <= number <= 17:
                    mouse = (number - 10) // 2 + 1
                    lesion = 1 if (number - 10) % 2 == 0 else 2
                    assignment = (
                        f"UVR_{mouse}", "UVR", mouse + 4,
                        f"{mouse}-AK-{lesion}",
                    )
                    role = f"AK-{lesion}"
            elif cscc_match:
                number = int(cscc_match.group(1))
                if 18 <= number <= 21:
                    mouse = number - 17
                    assignment = (
                        f"UVR_{mouse}", "UVR", mouse + 4, f"{mouse}-cSCC",
                    )
                    role = "cSCC"
            if not assignment:
                continue
        animal, exposure, subject, paper_name = assignment
        result.at[index, "publication_candidate_animal_id"] = animal
        result.at[index, "publication_candidate_exposure_group"] = exposure
        result.at[index, "publication_candidate_subject"] = subject
        result.at[index, "publication_candidate_sample_name"] = paper_name
        result.at[index, "publication_candidate_sample_role"] = role
        result.at[index, "publication_crosswalk_status"] = "provisional_publication_order_match"

    candidate = result.loc[
        result["publication_candidate_animal_id"].notna()
    ].copy()
    design_rows = []
    exposed = candidate.loc[
        candidate["publication_candidate_exposure_group"].eq("UVR")
    ]
    for animal, frame in exposed.groupby("publication_candidate_animal_id"):
        n_normal = int((frame["stage"] == "Normal skin").sum())
        n_ak = int((frame["stage"] == "Actinic keratosis").sum())
        n_cscc = int((frame["stage"] == "Cutaneous SCC").sum())
        design_rows.append({
            "candidate_animal_id": animal,
            "candidate_subject": int(frame["publication_candidate_subject"].iloc[0]),
            "n_normal_uvr": n_normal,
            "n_ak_lesions": n_ak,
            "n_cscc": n_cscc,
            "complete_candidate_matched_set": (
                n_normal == 1 and n_ak == 2 and n_cscc == 1
            ),
            "verification_status": "provisional_not_for_inference",
        })
    return result, pd.DataFrame(design_rows, columns=[
        "candidate_animal_id", "candidate_subject", "n_normal_uvr",
        "n_ak_lesions", "n_cscc", "complete_candidate_matched_set",
        "verification_status",
    ])


def experimental_design(metadata: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Assess whether animal-level repeated measures can be reconstructed safely."""
    animal_rows = metadata.loc[metadata["animal_id"].notna()].copy()
    if animal_rows.empty:
        return pd.DataFrame(columns=[
            "animal_id", "n_normal_uvr", "n_ak_lesions", "n_cscc", "complete_matched_set"
        ]), {
            "animal_id_available": False,
            "n_animals_identified": 0,
            "design_ready_for_paired_model": False,
            "inference_status": "blocked_missing_animal_id",
            "design_notice": (
                "Animal identifiers were not explicitly recoverable from GEO sample "
                "characteristics. Sample title suffixes are not treated as animal IDs."
            ),
        }
    exposed = animal_rows.loc[
        animal_rows["exposure_group"].eq("UVR")
        | animal_rows["stage"].isin(["Actinic keratosis", "Cutaneous SCC"])
    ]
    rows = []
    for animal, frame in exposed.groupby("animal_id"):
        n_normal = int((frame["analysis_group"] == "Normal skin (UVR)").sum())
        n_ak = int((frame["stage"] == "Actinic keratosis").sum())
        n_cscc = int((frame["stage"] == "Cutaneous SCC").sum())
        rows.append({
            "animal_id": animal,
            "n_normal_uvr": n_normal,
            "n_ak_lesions": n_ak,
            "n_cscc": n_cscc,
            "complete_matched_set": n_normal >= 1 and n_ak >= 1 and n_cscc >= 1,
        })
    design = pd.DataFrame(rows)
    complete = bool(not design.empty and design["complete_matched_set"].all())
    return design, {
        "animal_id_available": True,
        "n_animals_identified": int(design["animal_id"].nunique()),
        "design_ready_for_paired_model": complete,
        "inference_status": "design_reconstructed_model_not_yet_fitted" if complete
        else "blocked_incomplete_animal_stage_mapping",
        "design_notice": (
            "Animal-level links were extracted from explicit sample characteristics. "
            "Repeated AK lesions must be handled within animal; no differential-expression "
            "model is fitted by this workflow."
        ),
    }


def summarize_geo_counts(
    counts: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return library QC and descriptive target expression without inferential claims."""
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("At least one GEO library has zero total counts.")
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 0.5)
    qc = metadata.copy()
    qc["library_size"] = qc["sample_accession"].map(library_sizes).astype(int)
    qc["detected_features"] = qc["sample_accession"].map((counts > 0).sum()).astype(int)
    group_medians = qc.groupby("analysis_group")["library_size"].transform("median")
    qc["library_size_ratio_to_group_median"] = qc["library_size"] / group_medians
    qc["qc_flag"] = np.where(
        qc["library_size_ratio_to_group_median"] < 0.5,
        "review_low_library_size", "pass",
    )

    canonical = {str(gene).casefold(): gene for gene in TARGETS}
    matches = {
        feature: canonical[feature.casefold()]
        for feature in log_cpm.index if feature.casefold() in canonical
    }
    rows = []
    for feature, gene in matches.items():
        for group, sample_rows in metadata.groupby("analysis_group"):
            samples = sample_rows["sample_accession"].tolist()
            values = log_cpm.loc[feature, samples].astype(float)
            rows.append({
                "Gene": gene,
                "Module": TARGETS[gene],
                "Analysis_group": group,
                "Stage": sample_rows["stage"].iloc[0],
                "n_samples": len(samples),
                "median_log2_cpm": float(values.median()),
                "q1_log2_cpm": float(values.quantile(.25)),
                "q3_log2_cpm": float(values.quantile(.75)),
            })
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    contrasts = []
    if not summary.empty:
        for gene, frame in summary.groupby("Gene"):
            values = frame.set_index("Analysis_group")["median_log2_cpm"]
            for left, right in (
                ("Normal skin (UVF)", "Normal skin (UVR)"),
                ("Normal skin (UVR)", "Actinic keratosis"),
                ("Actinic keratosis", "Cutaneous SCC"),
                ("Normal skin (UVF)", "Cutaneous SCC"),
            ):
                if left in values and right in values:
                    contrasts.append({
                        "Gene": gene,
                        "Module": frame["Module"].iloc[0],
                        "contrast": f"{right}_vs_{left}",
                        "median_log2_cpm_difference": float(values[right] - values[left]),
                        "interpretation": (
                            "Descriptive difference of group medians; not differential-expression inference"
                        ),
                    })
    return qc, summary, pd.DataFrame(contrasts, columns=CONTRAST_COLUMNS)


def run_transcriptomics(
    output_dir: Path,
    client: NCBIGEOClient | None = None,
    cache_dir: Path | None = None,
) -> dict[str, Any]:
    """Retrieve GSE158634 from NCBI, validate the sample design and export summaries."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = cache_dir or RAW_DIR
    status: dict[str, Any] = {
        "accession": "GSE158634",
        "organism": "Mus musculus",
        "assay": "bulk RNA-seq",
        "source": "NCBI GEO E-utilities + official GEO FTP",
        "status": "not_run",
        "n_samples_expected": len(EXPECTED_SAMPLES),
        "inference_status": "not_performed",
        "notice": (
            "Independent mouse cohort; it is never merged with the human mutation cohort. "
            "Descriptive log2-CPM is not differential-expression inference."
        ),
    }
    try:
        client = client or NCBIGEOClient()
        resources = _download_geo_sources(cache_dir, client)
        soft_text = read_gzip_text(resources["soft_path"])
        metadata, series = parse_geo_family_soft(soft_text)
        sra_ids = metadata["sra_experiment_accession"].dropna().astype(str).tolist()
        sra_error = None
        try:
            sra_records = client.fetch_sra_experiments(sra_ids)
        except NCBIError as exc:
            LOGGER.warning("Linked SRA enrichment unavailable: %s", exc)
            sra_records = {}
            sra_error = str(exc)
        metadata = merge_sra_metadata(metadata, sra_records)
        biosample_ids = metadata["biosample_accession"].dropna().astype(str).tolist()
        biosample_error = None
        try:
            biosample_records = client.fetch_biosample_records(biosample_ids)
        except NCBIError as exc:
            LOGGER.warning("Linked BioSample enrichment unavailable: %s", exc)
            biosample_records = {}
            biosample_error = str(exc)
        metadata = merge_biosample_metadata(metadata, biosample_records)
        metadata, candidate_design = build_publication_candidate_crosswalk(metadata)
        observed = set(metadata["sample_accession"])
        missing = sorted(EXPECTED_SAMPLES - observed)
        unexpected = sorted(observed - EXPECTED_SAMPLES)
        if missing or unexpected:
            raise ValueError(
                f"GEO SOFT sample set mismatch; missing={missing}, unexpected={unexpected}"
            )
        counts, metadata = load_geo_counts(resources["counts_path"], metadata)
        if set(counts.columns) != EXPECTED_SAMPLES:
            raise ValueError("GEO count matrix does not contain the complete expected sample set.")
        qc, summary, contrasts = summarize_geo_counts(counts, metadata)
        design, design_status = experimental_design(metadata)

        metadata.to_csv(output_dir / "geo_sample_metadata.csv", index=False)
        design.to_csv(output_dir / "geo_experimental_design.csv", index=False)
        metadata.loc[metadata["publication_candidate_animal_id"].notna()].to_csv(
            output_dir / "geo_publication_candidate_crosswalk.csv", index=False
        )
        candidate_design.to_csv(output_dir / "geo_publication_candidate_design.csv", index=False)
        qc.to_csv(output_dir / "geo_sample_qc.csv", index=False)
        summary.to_csv(output_dir / "geo_target_expression.csv", index=False)
        contrasts.to_csv(output_dir / "geo_target_contrasts.csv", index=False)
        status.update({
            "status": "source_data_validated",
            "n_samples": int(counts.shape[1]),
            "n_features": int(counts.shape[0]),
            "n_target_genes": int(summary["Gene"].nunique()) if not summary.empty else 0,
            "n_qc_flagged": int(qc["qc_flag"].ne("pass").sum()),
            "sample_metadata_source": "GSE158634_family.soft.gz",
            "count_matrix_source": "GSE158634_3.count.csv.gz",
            "series_title": series.get("title"),
            "series_overall_design": " ".join(series.get("overall_design", [])),
            "sample_set_complete": True,
            "counts_are_nonnegative_integers": True,
            "feature_identifier_type": (
                "gene_symbol" if summary["Gene"].nunique() else "unresolved"
            ),
            **_metadata_consistency(metadata),
            "n_biosample_records": len(biosample_records),
            "biosample_enrichment_error": biosample_error,
            "n_sra_experiments_retrieved": len(sra_records),
            "sra_enrichment_error": sra_error,
            "sra_instrument_models": sorted(set(metadata["sra_instrument_model"].dropna().astype(str))),
            "normal_exposure_unknown_samples": int((metadata["analysis_group"] == "Normal skin (exposure unknown)").sum()),
            "normal_samples_split_by_exposure": bool(metadata.loc[metadata["stage"].eq("Normal skin"), "exposure_group"].notna().all()),
            "analysis_group_counts": metadata["analysis_group"].value_counts().to_dict(),
            "publication_design_reference": "PMID:33664254; eight mice split into UVF/UVR groups, with repeated tissues collected from UVR mice",
            "publication_candidate_crosswalk_samples": int(metadata["publication_candidate_animal_id"].notna().sum()),
            "publication_candidate_crosswalk_verified": False,
            "publication_candidate_uvr_animals_complete": bool(not candidate_design.empty and candidate_design["complete_candidate_matched_set"].all()),
            "publication_candidate_crosswalk_notice": "Provisional alignment of GEO title numbering to publication Table 2; do not use as verified animal metadata or for differential-expression inference.",
            "entrez_uid_count": len(resources.get("entrez_uids", [])),
            **design_status,
            "inference_status": design_status["inference_status"],
            "notice": (
                "NCBI sample metadata and integer count matrix were retrieved and reconciled. "
                "The deposited SOFT records do not currently expose animal IDs or UVF/UVR labels "
                "for each sample. The publication describes eight mice split between UVF and UVR, "
                "but sample-to-animal links are not accepted from title order as verified metadata. "
                "A separate provisional crosswalk is exported for review only. Descriptive expression "
                "only; no DE model is run until the author-confirmed animal-level design is reconstructed."
            ),
        })
    except (NCBIError, OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        LOGGER.exception("GSE158634 NCBI integration failed")
        status.update({
            "status": "retrieval_or_validation_failed",
            "validation_message": str(exc),
            "inference_status": "not_performed",
        })
        for filename, columns in (
            ("geo_sample_metadata.csv", ["sample_accession", "sample_title"]),
            ("geo_experimental_design.csv", ["animal_id", "complete_matched_set"]),
            ("geo_publication_candidate_crosswalk.csv", ["sample_accession", "sample_title", "publication_candidate_animal_id", "publication_crosswalk_status"]),
            ("geo_publication_candidate_design.csv", ["candidate_animal_id", "complete_candidate_matched_set"]),
            ("geo_sample_qc.csv", ["sample_accession", "stage", "library_size", "qc_flag"]),
            ("geo_target_expression.csv", SUMMARY_COLUMNS),
            ("geo_target_contrasts.csv", CONTRAST_COLUMNS),
        ):
            pd.DataFrame(columns=columns).to_csv(output_dir / filename, index=False)
    (output_dir / "geo_transcriptomics_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return status
