"""Prepare aggregate-only JSON for the public GitHub Pages dashboard."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from akstats.analysis import gene_level_maf

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs"
DESTINATION = ROOT / "docs" / "data" / "dashboard.json"


def records(filename: str) -> list[dict]:
    path = OUTPUT / filename
    if not path.exists():
        raise FileNotFoundError(f"Required research output is missing: {path}")
    frame = pd.read_csv(path)
    return json.loads(frame.to_json(orient="records"))


def main() -> None:
    overview = pd.read_excel(OUTPUT / "AKStats_research_results.xlsx", sheet_name="Study overview")
    metrics = dict(zip(overview["Metric"], overview["Value"]))
    maf = records("maf_descriptive.csv")
    modules = records("module_summary.csv")
    fisher = records("fisher_tests.csv")
    patient_metrics = records("patient_level_summary.csv")
    patient_tests = records("patient_level_tests.csv")
    geo_status = json.loads((OUTPUT / "geo_transcriptomics_status.json").read_text(encoding="utf-8"))
    geo_qc = records("geo_sample_qc.csv")
    geo_expression = records("geo_target_expression.csv")
    geo_contrasts = records("geo_target_contrasts.csv")

    # Keep each analytical level distinct: group-level MAF, module-level ΔMAF,
    # and gene-level MAF. Never assign one arbitrary module to an entire group.
    axis = pd.read_csv(ROOT / "data" / "processed" / "arid2_mapk_variants.csv")
    gene_maf = gene_level_maf(axis)

    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "panel_genes": int(metrics["Panel genes"]),
            "source_variants": int(metrics["Source variant records"]),
            "axis_variants": int(metrics["Filtered functional axis variants"]),
            "patients": int(metrics["Patients in source table"]),
        },
        "maf": maf,
        "modules": modules,
        "gene_maf": json.loads(gene_maf.to_json(orient="records")),
        "fisher": fisher,
        "patient_metrics": patient_metrics,
        "patient_tests": patient_tests,
        "geo_transcriptomics": geo_status,
        "geo_sample_qc": geo_qc,
        "geo_expression": geo_expression,
        "geo_contrasts": geo_contrasts,
        "data_notice": "Aggregate exploratory molecular summaries; no individual clinical records or diagnostic outputs.",
    }
    DESTINATION.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
