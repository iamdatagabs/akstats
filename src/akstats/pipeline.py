"""End-to-end mutation analysis orchestration."""
import json
import logging

from akstats.analysis import (
    cooccurrence_table,
    descriptive_maf,
    fisher_by_stage,
    module_summary,
    paired_stage_tests,
    prepare_axis,
)
from akstats.config import OUTPUT_DIR, PROCESSED_DIR
from akstats.data import load_mutation_sources
from akstats.patient_metrics import patient_level_maf, patient_level_summary, patient_level_tests
from akstats.plots import save_figures
from akstats.report import export_research_products

LOGGER = logging.getLogger(__name__)


def run_pipeline() -> None:
    panel, variants = load_mutation_sources()
    axis = prepare_axis(panel, variants)
    if axis.empty:
        raise RuntimeError("No functional ARID2/MAPK variants matched the sequenced panel.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    axis.to_csv(PROCESSED_DIR / "arid2_mapk_variants.csv", index=False)
    descriptive_maf(axis).to_csv(OUTPUT_DIR / "maf_descriptive.csv", index=False)
    module_summary(axis).to_csv(OUTPUT_DIR / "module_summary.csv", index=False)
    cooc = cooccurrence_table(axis, variants["Case"])
    cooc.to_csv(OUTPUT_DIR / "patient_cooccurrence.csv", index=False)
    fisher_by_stage(cooc).to_csv(OUTPUT_DIR / "fisher_tests.csv", index=False)
    patient_maf = patient_level_maf(axis, variants["Case"])
    patient_level_summary(patient_maf).to_csv(OUTPUT_DIR / "patient_level_summary.csv", index=False)
    patient_level_tests(patient_maf).to_csv(OUTPUT_DIR / "patient_level_tests.csv", index=False)
    tests = paired_stage_tests(axis)
    (OUTPUT_DIR / "exploratory_tests.json").write_text(
        json.dumps(tests, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    figures = save_figures(axis, OUTPUT_DIR)
    products = export_research_products(OUTPUT_DIR, panel, variants, axis, figures)
    LOGGER.info("Panel genes: %d", panel["Gene"].nunique())
    LOGGER.info(
        "Raw variants: %d | filtered axis variants: %d | patients: %d",
        len(variants), len(axis), variants["Case"].nunique(),
    )
    LOGGER.info("Research products: %s", [str(path) for path in products])
