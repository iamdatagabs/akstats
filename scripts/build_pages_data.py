"""Prepare the public, aggregate-only JSON consumed by the GitHub Pages dashboard."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

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
    # Module labels are assigned from the variant-level descriptive summaries.
    axis = pd.read_csv(ROOT / "data" / "processed" / "arid2_mapk_variants.csv")
    module_map = axis[["Biological_group", "Module"]].drop_duplicates()
    for row in maf:
        candidates = module_map.loc[module_map["Biological_group"].eq(row["Biological_group"]), "Module"]
        row["Module"] = candidates.iloc[0] if not candidates.empty else row["Biological_group"]
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
        "fisher": fisher,
        "patient_metrics": patient_metrics,
        "patient_tests": patient_tests,
        "data_notice": "Aggregate exploratory molecular summaries; no individual clinical records or diagnostic outputs.",
    }
    DESTINATION.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
