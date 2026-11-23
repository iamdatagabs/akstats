"""Publication-oriented research product exporters."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import markdown
import pandas as pd


def export_research_products(output_dir: Path, panel: pd.DataFrame, variants: pd.DataFrame,
                             axis: pd.DataFrame, figures: list[Path]) -> list[Path]:
    """Create a consolidated workbook, Markdown report, HTML report and manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "Filtered variants": axis,
        "MAF descriptive": pd.read_csv(output_dir / "maf_descriptive.csv"),
        "Pathway modules": pd.read_csv(output_dir / "module_summary.csv"),
        "Patient cooccurrence": pd.read_csv(output_dir / "patient_cooccurrence.csv"),
        "Fisher tests": pd.read_csv(output_dir / "fisher_tests.csv"),
        "Patient-level summary": pd.read_csv(output_dir / "patient_level_summary.csv"),
        "Patient-level tests": pd.read_csv(output_dir / "patient_level_tests.csv"),
        "GEO sample metadata": pd.read_csv(output_dir / "geo_sample_metadata.csv"),
        "GEO experimental design": pd.read_csv(output_dir / "geo_experimental_design.csv"),
        "GEO candidate crosswalk": pd.read_csv(output_dir / "geo_publication_candidate_crosswalk.csv"),
        "GEO candidate design": pd.read_csv(output_dir / "geo_publication_candidate_design.csv"),
        "GEO sample QC": pd.read_csv(output_dir / "geo_sample_qc.csv"),
        "GEO target expression": pd.read_csv(output_dir / "geo_target_expression.csv"),
        "GEO target contrasts": pd.read_csv(output_dir / "geo_target_contrasts.csv"),
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
    patient_summary, patient_tests = tables["Patient-level summary"], tables["Patient-level tests"]
    geo_status_path = output_dir / "geo_transcriptomics_status.json"
    geo_status = json.loads(geo_status_path.read_text(encoding="utf-8")) if geo_status_path.exists() else {}
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
        "## Patient-level molecular summary", "", patient_summary.to_markdown(index=False), "",
        "The patient-level MAF measure is the maximum variant MAF observed within each patient, group and stage. Patients without an eligible variant are represented as zero in the filtered dataset; this does not mean a biological MAF measurement of zero.", "",
        "## Patient-level paired tests", "", patient_tests.to_markdown(index=False), "",
        "Friedman compares the three paired tissue stages; Wilcoxon compares AK with cSCC; exact McNemar compares threshold-positive status. These analyses are exploratory and do not establish clinical utility.", "",
        "## Patient-level co-occurrence", "",
        f"Patients in co-occurrence matrix: **{len(cooc)}**", "", cooc.to_markdown(index=False), "",
        "## Fisher exact tests", "", fisher.to_markdown(index=False), "",
        "## Figures", "",
    ]
    lines.extend(f"- [{path.name}]({path.name})" for path in figures)
    lines.extend([
        "", "## Independent mouse transcriptomics — GSE158634", "",
        f"- NCBI retrieval/validation status: **{geo_status.get('status', 'not_run')}**",
        f"- Samples reconciled: **{geo_status.get('n_samples', 0)} / 20**",
        f"- Target gene symbols mapped: **{geo_status.get('n_target_genes', 0)}**",
        f"- Libraries flagged for review: **{geo_status.get('n_qc_flagged', 0)}**",
        f"- Animal IDs available: **{geo_status.get('animal_id_available', 'unknown')}**",
        f"- Experimental design: **{geo_status.get('inference_status', 'not assessed')}**",
        f"- Genome-build consistency: **{geo_status.get('genome_build_consistency', 'not assessed')}**",
        f"- Genome builds mentioned: **{', '.join(geo_status.get('genome_builds_detected', [])) or 'none detected'}**",
        f"- Linked SRA experiments retrieved: **{geo_status.get('n_sra_experiments_retrieved', 0)}**",
        f"- SRA sample accessions reconciled: **{geo_status.get('n_sra_sample_accessions', 0)}**",
        f"- Linked BioSamples retrieved: **{geo_status.get('n_biosample_records', 0)}**",
        f"- Exposure metadata conflicts: **{geo_status.get('n_exposure_metadata_conflicts', 0)}**",
        f"- Animal-ID metadata conflicts: **{geo_status.get('n_animal_id_metadata_conflicts', 0)}**",
        f"- SRA instrument models: **{', '.join(geo_status.get('sra_instrument_models', [])) or 'not available'}**", "",
        f"- Instrument reported by article: **{geo_status.get('publication_instrument_reported', 'not reported')}**",
        geo_status.get("platform_instrument_warning") or "",
        geo_status.get("notice", "No transcriptomic inference was performed."), "",
        "### GEO sample metadata", "",
        tables["GEO sample metadata"].to_markdown(index=False)
        if not tables["GEO sample metadata"].empty else "No GEO sample metadata available.", "",
        "### Animal-level design audit", "",
        tables["GEO experimental design"].to_markdown(index=False)
        if not tables["GEO experimental design"].empty else "Verified animal-level links were not reconstructed.", "",
        "### Provisional publication-based crosswalk", "",
        tables["GEO candidate design"].to_markdown(index=False)
        if not tables["GEO candidate design"].empty else "No candidate animal groups could be aligned.", "",
        "The candidate crosswalk aligns GEO sample-title numbering with the order of sample names in Table 2 of PMID 33664254. It is a review aid, not an author-confirmed GSM-to-animal key. It is not copied into verified animal IDs and is not used for inferential modeling.", "",
        tables["GEO candidate crosswalk"].to_markdown(index=False)
        if not tables["GEO candidate crosswalk"].empty else "No provisional sample crosswalk available.", "",
        "### Sample-level QC", "",
        tables["GEO sample QC"].to_markdown(index=False)
        if not tables["GEO sample QC"].empty else "No GEO sample QC available.", "",
        "### Descriptive target expression", "",
        tables["GEO target expression"].to_markdown(index=False)
        if not tables["GEO target expression"].empty else "No target-gene expression table available.", "",
        "Contrasts are descriptive differences between group medians of log2(CPM + 0.5), not differential-expression tests. No sample is excluded solely by the library-size flag. No inferential p-values or enrichment results are reported until the animal-level experimental unit and model design are verified.", "",
        "## Interpretation and limitations", "",
        "- Exploratory re-analysis of published supplementary data; not a diagnostic or prognostic model.",
        "- Variant-level MAF observations are clustered within patients; row-wise Friedman/Wilcoxon p-values are exploratory, not patient-level inference.",
        "- ΔMAF is descriptive, not clone size. Purity, copy number, sampling and assay properties can influence MAF.",
        "- The 2% MAF threshold is operational, not a universal biological cutoff.",
        "- Co-occurrence does not establish biological interaction or causation. Interpret Fisher tests with cell counts and cohort size.",
        "- The mouse RNA-seq cohort is independent from the human mutation cohort and must never be joined by patient identifiers.",
        "- The GEO workflow retrieves sample characteristics from the official family SOFT file. Animal IDs are accepted only when explicitly present in deposited metadata; sample-title suffixes are not treated as verified animal IDs. A separate provisional publication crosswalk is provided for manual review only.",
        "- Genome assembly and annotation release must be reconciled. GSE158634 sample processing metadata currently mentions both hg19 and mm9; do not select one silently.",
        "- The GEO platform GPL18480 is titled HiSeq 1500, while the sample processing text mentions HiSeq 2500; verify the instrument against the publication.",
        "- Review the source article and annotation protocol before biological interpretation or publication.", "",
        "## Data provenance", "",
        "Mutation panel and variant tables: supplementary files MOESM8/MOESM9 associated with DOI 10.1038/s41467-025-65687-y. Original workbooks are versioned under data/raw/; processed outputs are generated separately.",
        "Transcriptomics metadata: NCBI GEO GSE158634 family SOFT via official GEO FTP; Entrez GEO DataSets E-utilities are queried for accession discovery/summary. Counts: GSE158634_3.count.csv.gz from the official GEO supplementary directory.", "",
        "NCBI E-utilities usage is subject to NCBI rate limits. Configure NCBI_EMAIL and optionally NCBI_API_KEY as environment variables; no API key is required for standard low-volume use.", "",
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
        "transcriptomics_source": "NCBI GEO GSE158634",
        "transcriptomics_status": geo_status.get("status", "not_run"),
        "panel_gene_count": panel["Gene"].nunique(),
        "source_variant_count": len(variants),
        "filtered_variant_count": len(axis),
        "patient_count": variants["Case"].nunique(),
        "processed_dataset": "data/processed/arid2_mapk_variants.csv",
        "figures": [path.name for path in figures],
        "outputs": [workbook.name, md_path.name, html_path.name],
    }, indent=2), encoding="utf-8")
    return [workbook, md_path, html_path, manifest]
