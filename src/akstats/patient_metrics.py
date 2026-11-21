"""Patient-level molecular metrics; patients, not variant rows, are the unit of analysis."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

from akstats.config import CLONALITY_THRESHOLD

STAGES = {"Normal": "Normal_MAF", "AK": "AK_MAF", "cSCC": "SCC_MAF"}


def patient_level_maf(axis: pd.DataFrame, all_cases: pd.Series) -> pd.DataFrame:
    """Summarize the maximum observed variant MAF per patient, group and stage.

    Zero means no qualifying variant was observed in the filtered axis for that
    patient/group/stage; it is not a measured wild-type MAF.
    """
    cases = sorted(all_cases.dropna().astype(str).unique())
    rows: list[dict] = []
    for case in cases:
        patient = axis.loc[axis["Case"].astype(str).eq(case)]
        for group in ("ARID2", "MAPK"):
            subset = patient.loc[patient["Biological_group"].eq(group)]
            for stage, column in STAGES.items():
                values = subset[column].dropna()
                rows.append({
                    "Case": case,
                    "Biological_group": group,
                    "Stage": stage,
                    "max_variant_maf": float(values.max()) if not values.empty else 0.0,
                    "variant_count": int(values.size),
                    "present_at_threshold": bool(values.ge(CLONALITY_THRESHOLD).any()),
                })
    return pd.DataFrame(rows)


def patient_level_summary(patient_maf: pd.DataFrame) -> pd.DataFrame:
    """Aggregate patient-level MAF and prevalence with Wilson 95% confidence intervals."""
    rows: list[dict] = []
    for (group, stage), frame in patient_maf.groupby(
        ["Biological_group", "Stage"], sort=False
    ):
        n = int(frame["Case"].nunique())
        positive = int(frame.drop_duplicates("Case")["present_at_threshold"].sum())
        low, high = proportion_confint(positive, n, method="wilson") if n else (np.nan, np.nan)
        values = frame.drop_duplicates("Case")["max_variant_maf"]
        rows.append({
            "Biological_group": group,
            "Stage": stage,
            "n_patients": n,
            "median_patient_max_maf": float(values.median()) if n else np.nan,
            "q1_patient_max_maf": float(values.quantile(.25)) if n else np.nan,
            "q3_patient_max_maf": float(values.quantile(.75)) if n else np.nan,
            "patients_positive": positive,
            "prevalence": positive / n if n else np.nan,
            "prevalence_ci_low": float(low),
            "prevalence_ci_high": float(high),
        })
    return pd.DataFrame(rows)


def patient_level_tests(patient_maf: pd.DataFrame) -> pd.DataFrame:
    """Paired stage comparisons using one observation per patient and group."""
    rows: list[dict] = []
    for group in ("ARID2", "MAPK"):
        frame = patient_maf.loc[patient_maf["Biological_group"].eq(group)]
        wide = frame.pivot(index="Case", columns="Stage", values="max_variant_maf")
        wide = wide.reindex(columns=["Normal", "AK", "cSCC"]).fillna(0.0)
        if len(wide) < 2:
            friedman_stat = friedman_p = wilcoxon_stat = wilcoxon_p = np.nan
        else:
            friedman = stats.friedmanchisquare(wide["Normal"], wide["AK"], wide["cSCC"])
            friedman_stat, friedman_p = float(friedman.statistic), float(friedman.pvalue)
            try:
                paired = stats.wilcoxon(wide["AK"], wide["cSCC"], zero_method="wilcox")
                wilcoxon_stat, wilcoxon_p = float(paired.statistic), float(paired.pvalue)
            except ValueError:
                wilcoxon_stat = wilcoxon_p = np.nan
        ak = frame.loc[frame["Stage"].eq("AK")].set_index("Case")["present_at_threshold"]
        scc = frame.loc[frame["Stage"].eq("cSCC")].set_index("Case")["present_at_threshold"]
        paired_presence = pd.concat([ak, scc], axis=1).fillna(False)
        gained = int((~paired_presence.iloc[:, 0] & paired_presence.iloc[:, 1]).sum())
        lost = int((paired_presence.iloc[:, 0] & ~paired_presence.iloc[:, 1]).sum())
        discordant = gained + lost
        mcnemar_p = float(stats.binomtest(min(gained, lost), discordant, .5).pvalue) if discordant else 1.0
        rows.append({
            "Biological_group": group,
            "n_patients": len(wide),
            "friedman_statistic": friedman_stat,
            "friedman_p": friedman_p,
            "wilcoxon_statistic": wilcoxon_stat,
            "wilcoxon_p": wilcoxon_p,
            "presence_gained_AK_to_cSCC": gained,
            "presence_lost_AK_to_cSCC": lost,
            "mcnemar_exact_p": mcnemar_p,
        })
    return pd.DataFrame(rows)
