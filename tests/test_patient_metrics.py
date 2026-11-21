import pandas as pd

from akstats.patient_metrics import (
    patient_level_maf,
    patient_level_summary,
    patient_level_tests,
)


def test_patient_level_metrics_use_patients_not_variant_rows():
    axis = pd.DataFrame([
        {"Case": "P1", "Biological_group": "ARID2", "Normal_MAF": 0.0, "AK_MAF": 0.03, "SCC_MAF": 0.12},
        {"Case": "P1", "Biological_group": "ARID2", "Normal_MAF": 0.0, "AK_MAF": 0.05, "SCC_MAF": 0.10},
        {"Case": "P1", "Biological_group": "MAPK", "Normal_MAF": 0.0, "AK_MAF": 0.01, "SCC_MAF": 0.08},
        {"Case": "P2", "Biological_group": "MAPK", "Normal_MAF": 0.0, "AK_MAF": 0.00, "SCC_MAF": 0.00},
        {"Case": "P3", "Biological_group": "ARID2", "Normal_MAF": 0.0, "AK_MAF": 0.00, "SCC_MAF": 0.04},
    ])
    cases = pd.Series(["P1", "P2", "P3"])
    patient = patient_level_maf(axis, cases)
    # One row per patient x group x stage; repeated variants do not inflate n.
    assert len(patient) == 3 * 2 * 3
    p1_ak = patient.query("Case == 'P1' and Biological_group == 'ARID2' and Stage == 'AK'")
    assert p1_ak["max_variant_maf"].iloc[0] == 0.05
    summary = patient_level_summary(patient)
    arid_ak = summary.query("Biological_group == 'ARID2' and Stage == 'AK'").iloc[0]
    assert arid_ak["n_patients"] == 3
    assert arid_ak["patients_positive"] == 1
    assert 0 <= arid_ak["prevalence_ci_low"] <= arid_ak["prevalence_ci_high"] <= 1


def test_patient_level_tests_return_one_row_per_group():
    axis = pd.DataFrame([
        {"Case": "P1", "Biological_group": "ARID2", "Normal_MAF": 0.0, "AK_MAF": 0.03, "SCC_MAF": 0.12},
        {"Case": "P2", "Biological_group": "ARID2", "Normal_MAF": 0.0, "AK_MAF": 0.00, "SCC_MAF": 0.04},
        {"Case": "P1", "Biological_group": "MAPK", "Normal_MAF": 0.0, "AK_MAF": 0.01, "SCC_MAF": 0.08},
        {"Case": "P2", "Biological_group": "MAPK", "Normal_MAF": 0.0, "AK_MAF": 0.00, "SCC_MAF": 0.00},
    ])
    patient = patient_level_maf(axis, pd.Series(["P1", "P2"]))
    tests = patient_level_tests(patient)
    assert set(tests["Biological_group"]) == {"ARID2", "MAPK"}
    assert (tests["n_patients"] == 2).all()
    assert tests["mcnemar_exact_p"].between(0, 1).all()
