"""Transformations and exploratory statistical summaries."""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from akstats.config import CLONALITY_THRESHOLD, FUNCTIONAL_CLASSES, GENE_TO_MODULE, MAF_COLUMNS


def prepare_axis(panel: pd.DataFrame, variants: pd.DataFrame) -> pd.DataFrame:
    """Annotate variants and retain functional ARID2/MAPK genes in the sequenced panel."""
    sequenced = set(panel["Gene"].dropna().astype(str))
    allowed = {gene for gene in GENE_TO_MODULE if gene in sequenced}
    data = variants.copy()
    data["Hugo_Symbol"] = data["Hugo_Symbol"].astype(str)
    data["Module"] = data["Hugo_Symbol"].map(GENE_TO_MODULE)
    data["Biological_group"] = np.where(data["Hugo_Symbol"].eq("ARID2"), "ARID2", "MAPK")
    data["Is_functional"] = data["Variant_Classification"].isin(FUNCTIONAL_CLASSES)
    data["Delta_MAF"] = data["SCC_MAF"] - data["AK_MAF"]
    data["Clonal_dynamics"] = np.select(
        [
            data["AK_MAF"].lt(CLONALITY_THRESHOLD) & data["SCC_MAF"].ge(CLONALITY_THRESHOLD),
            data["AK_MAF"].ge(CLONALITY_THRESHOLD) & data["SCC_MAF"].ge(CLONALITY_THRESHOLD),
            data["AK_MAF"].ge(CLONALITY_THRESHOLD) & data["SCC_MAF"].lt(CLONALITY_THRESHOLD),
        ],
        ["Emergent at cSCC", "Detected at both stages", "Not detected above threshold at cSCC"],
        default="Below threshold at both stages",
    )
    return data.loc[data["Hugo_Symbol"].isin(allowed) & data["Is_functional"]].copy()


def descriptive_maf(axis: pd.DataFrame) -> pd.DataFrame:
    """Summarize MAF by biological group and tissue stage."""
    rows = []
    for group, frame in axis.groupby("Biological_group", observed=True):
        for column in MAF_COLUMNS:
            values = frame[column].dropna()
            rows.append({
                "Biological_group": group, "Stage": column.removesuffix("_MAF"),
                "n_variants": int(values.size),
                "mean": float(values.mean()) if len(values) else np.nan,
                "sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
                "median": float(values.median()) if len(values) else np.nan,
                "q1": float(values.quantile(.25)) if len(values) else np.nan,
                "q3": float(values.quantile(.75)) if len(values) else np.nan,
            })
    return pd.DataFrame(rows)


def module_summary(axis: pd.DataFrame) -> pd.DataFrame:
    """Summarize variant-level MAF changes across pathway modules."""
    if axis.empty:
        return pd.DataFrame()
    return axis.groupby("Module", observed=True).agg(
        n_variants=("Hugo_Symbol", "size"),
        genes=("Hugo_Symbol", lambda x: ", ".join(sorted(x.unique()))),
        median_ak_maf=("AK_MAF", "median"),
        median_scc_maf=("SCC_MAF", "median"),
        median_delta_maf=("Delta_MAF", "median"),
        max_scc_maf=("SCC_MAF", "max"),
    ).reset_index()


def paired_stage_tests(axis: pd.DataFrame) -> dict[str, float]:
    """Run row-wise tests retained from the notebook; results are exploratory only.

    Rows are variants, not independent patients. Repeated variants per patient violate
    the independence assumption, so these p-values are not patient-level inference.
    """
    if len(axis) < 2:
        return {"friedman_statistic": np.nan, "friedman_p": np.nan,
                "wilcoxon_statistic": np.nan, "wilcoxon_p": np.nan}
    friedman = stats.friedmanchisquare(axis["Normal_MAF"], axis["AK_MAF"], axis["SCC_MAF"])
    wilcoxon = stats.wilcoxon(axis["AK_MAF"], axis["SCC_MAF"], zero_method="wilcox")
    return {
        "friedman_statistic": float(friedman.statistic), "friedman_p": float(friedman.pvalue),
        "wilcoxon_statistic": float(wilcoxon.statistic), "wilcoxon_p": float(wilcoxon.pvalue),
    }


def cooccurrence_table(axis: pd.DataFrame, all_cases: pd.Series) -> pd.DataFrame:
    """Build patient-level ARID2/MAPK presence at the configured MAF threshold."""
    rows = []
    for case in sorted(all_cases.dropna().unique()):
        sample = axis.loc[axis["Case"].eq(case)]
        row = {"Case": case}
        for stage, column in (("AK", "AK_MAF"), ("cSCC", "SCC_MAF")):
            row[f"{stage}_ARID2"] = bool(
                (sample["Hugo_Symbol"].eq("ARID2") & sample[column].ge(CLONALITY_THRESHOLD)).any()
            )
            row[f"{stage}_MAPK"] = bool(
                (sample["Hugo_Symbol"].ne("ARID2") & sample[column].ge(CLONALITY_THRESHOLD)).any()
            )
        rows.append(row)
    return pd.DataFrame(rows)


def fisher_by_stage(cooccurrence: pd.DataFrame) -> pd.DataFrame:
    """Run two-sided Fisher tests with explicit 2x2 cell counts."""
    rows = []
    for stage in ("AK", "cSCC"):
        a = cooccurrence[f"{stage}_ARID2"].astype(bool)
        b = cooccurrence[f"{stage}_MAPK"].astype(bool)
        table = np.array([
            [int((a & b).sum()), int((a & ~b).sum())],
            [int((~a & b).sum()), int((~a & ~b).sum())],
        ])
        odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
        rows.append({
            "stage": stage, "both_present": table[0, 0], "arid2_only": table[0, 1],
            "mapk_only": table[1, 0], "neither": table[1, 1],
            "odds_ratio": odds_ratio, "p_value": p_value,
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["p_fdr_bh"] = multipletests(result["p_value"], method="fdr_bh")[1]
    return result
