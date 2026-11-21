"""Publication-oriented research product exporters."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import markdown
import pandas as pd


def export_research_products(output_dir: Path, panel: pd.DataFrame, variants: pd.DataFrame,
                             axis: pd.DataFrame, figures: list[Path]) -> list[Path]:
    """Create a consolidated XLSX workbook, Markdown report, HTML report and manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "Filtered variants": axis,
        "MAF descriptive": pd.read_csv(output_dir / "maf_descriptive.csv"),
        "Pathway modules": pd.read_csv(output_dir / "module_summary.csv"),
        "Patient cooccurrence": pd.read_csv(output_dir / "patient_cooccurrence.csv"),
        "Fisher tests": pd.read_csv(output_dir / "fisher_tests.csv"),
        "Patient-level summary": pd.read_csv(output_dir / "patient_level_summary.csv"),
        "Patient-level tests": pd.read_csv(output_dir / "patient_level_tests.csv"),
    }
    workbook = output_dir / "AKStats_research_results.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        for sheet, frame in tables.items():
            frame.to_excel(writer, sheet_name=sheet[:31], index=False)
        pd.DataFrame([
            {"Metric": "Panel genes", "Value": panel["Gene"].nunique()},
            {"Metric": "Source variant records", "Value": len(variants)},
            {"Metric": "Filtered functional axis variants", "Value": len(axis)},
            {"Metric": "Patients in source table", "Value": variants["Case"].nunique()},
        ]).to_excel(writer, sheet_name="Study overview", index=False)

    descriptive, modules = tables["MAF descriptive"], tables["Pathway modules"]
    cooc, fisher = tables["Patient cooccurrence"], tables["Fisher tests"]
    lines = [
        "# AKStats — research results", "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat(timespec='seconds')}", "",
        "## Dataset overview", "",
        f"- Genes in sequenced panel: **{panel['Gene'].nunique()}**",
        f"- Variant records in source table: **{len(variants)}**",
        f"- Functional ARID2/MAPK variants retained: **{len(axis)}**",
        f"- Patients represented in source table: **{variants['Case'].nunique()}**", "",
        "## Descriptive MAF summary", "", descriptive.to_markdown(index=False), "",
        "## Pathway-module summary", "",
        modules.to_markdown(index=False) if not modules.empty else "No module-level records.", "",
        "## Patient-level co-occurrence", "",
        f"Patients in co-occurrence matrix: **{len(cooc)}**", "", cooc.to_markdown(index=False), "",
        "## Fisher exact tests", "", fisher.to_markdown(index=False), "", "## Figures", "",
    ]
    lines.extend(f"- [{path.name}]({path.name})" for path in figures)
    lines.extend([
        "", "## Interpretation and limitations", "",
        "- Exploratory re-analysis of published supplementary data; not a diagnostic or prognostic model.",
        "- Variant-level MAF observations are clustered within patients; row-wise Friedman/Wilcoxon p-values are exploratory, not patient-level inference.",
        "- ΔMAF is descriptive, not clone size. Purity, copy number, sampling and assay properties can influence MAF.",
        "- The 2% MAF threshold is operational, not a universal biological cutoff.",
        "- Co-occurrence does not establish biological interaction or causation. Interpret Fisher tests with cell counts and cohort size.",
        "- Review the source article and annotation protocol before biological interpretation or publication.", "",
        "## Data provenance", "",
        "Mutation panel and variant tables: supplementary files MOESM8/MOESM9 associated with DOI 10.1038/s41467-025-65687-y. Original workbooks are versioned under data/raw/; processed outputs are generated separately.", "",
    ])
    markdown_text = "\n".join(lines)
    md_path = output_dir / "research_report.md"
    md_path.write_text(markdown_text, encoding="utf-8")
    html_path = output_dir / "research_report.html"
    html_body = markdown.markdown(markdown_text, extensions=["tables"])
    html_path.write_text(
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<title>AKStats research results</title><style>"
        "body{font:16px/1.6 system-ui,sans-serif;max-width:1100px;margin:3rem auto;padding:0 1.2rem;color:#29251f}"
        "h1,h2{color:#594536}table{border-collapse:collapse;width:100%;display:block;overflow-x:auto}"
        "th,td{border:1px solid #ddd;padding:.5rem;text-align:left}th{background:#eee8df}"
        "a{color:#526b48}img{max-width:100%;height:auto}</style></head><body>"
        + html_body + "</body></html>", encoding="utf-8"
    )
    manifest = output_dir / "run_manifest.json"
    manifest.write_text(json.dumps({
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_doi": "10.1038/s41467-025-65687-y",
        "source_files": [
            "data/raw/41467_2025_65687_MOESM8_ESM.xlsx",
            "data/raw/41467_2025_65687_MOESM9_ESM.xlsx",
        ],
        "panel_gene_count": panel["Gene"].nunique(),
        "source_variant_count": len(variants),
        "filtered_variant_count": len(axis),
        "patient_count": variants["Case"].nunique(),
        "processed_dataset": "data/processed/arid2_mapk_variants.csv",
        "figures": [path.name for path in figures],
        "outputs": [workbook.name, md_path.name, html_path.name],
    }, indent=2), encoding="utf-8")
    return [workbook, md_path, html_path, manifest]
