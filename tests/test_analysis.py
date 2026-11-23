import numpy as np
import pandas as pd

from akstats.analysis import (
    cooccurrence_table,
    descriptive_maf,
    fisher_by_stage,
    gene_level_maf,
    prepare_axis,
)


def fixtures():
    panel = pd.DataFrame({"Gene": ["ARID2", "HRAS", "EGFR", "TP53"]})
    variants = pd.DataFrame([
        {"Case": "P01", "Hugo_Symbol": "ARID2", "Variant_Classification": "Missense_Mutation",
         "Normal_MAF": 0.0, "AK_MAF": 0.03, "SCC_MAF": 0.12},
        {"Case": "P01", "Hugo_Symbol": "HRAS", "Variant_Classification": "Missense_Mutation",
         "Normal_MAF": 0.0, "AK_MAF": 0.01, "SCC_MAF": 0.22},
        {"Case": "P02", "Hugo_Symbol": "TP53", "Variant_Classification": "Missense_Mutation",
         "Normal_MAF": 0.0, "AK_MAF": 0.10, "SCC_MAF": 0.20},
        {"Case": "P02", "Hugo_Symbol": "EGFR", "Variant_Classification": "Silent",
         "Normal_MAF": 0.0, "AK_MAF": 0.0, "SCC_MAF": 0.0},
    ])
    return panel, variants


def test_axis_filtering_and_delta():
    panel, variants = fixtures()
    result = prepare_axis(panel, variants)
    assert set(result["Hugo_Symbol"]) == {"ARID2", "HRAS"}
    np.testing.assert_allclose(result["Delta_MAF"], [0.09, 0.21])


def test_descriptive_summary_contains_each_stage():
    panel, variants = fixtures()
    result = descriptive_maf(prepare_axis(panel, variants))
    assert len(result) == 6
    assert set(result["Stage"]) == {"Normal", "AK", "SCC"}


def test_cooccurrence_is_patient_level():
    panel, variants = fixtures()
    axis = prepare_axis(panel, variants)
    cooc = cooccurrence_table(axis, variants["Case"])
    assert len(cooc) == 2
    assert bool(cooc.loc[cooc["Case"].eq("P01"), "AK_ARID2"].iloc[0])
    assert len(fisher_by_stage(cooc)) == 2


def test_gene_level_maf_preserves_gene_and_module_levels():
    panel, variants = fixtures()
    axis = prepare_axis(panel, variants)
    result = gene_level_maf(axis)
    assert set(result["Gene"]) == {"ARID2", "HRAS"}
    assert set(result.loc[result["Gene"].eq("ARID2"), "Module"]) == {"0. Chromatin (ARID2)"}
    assert set(result.loc[result["Gene"].eq("HRAS"), "Module"]) == {"2. RAS GTPases"}
    assert set(result["Stage"]) == {"Normal", "AK", "SCC"}
    assert result.groupby(["Gene", "Stage"]).size().eq(1).all()
