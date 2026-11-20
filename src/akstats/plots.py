"""Exploratory figure generation."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def save_figures(axis: pd.DataFrame, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    paths = []
    if axis.empty:
        return paths
    long = axis.melt(
        id_vars=["Case", "Hugo_Symbol", "Biological_group"],
        value_vars=["Normal_MAF", "AK_MAF", "SCC_MAF"],
        var_name="Stage", value_name="MAF",
    )
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.boxplot(data=long, x="Stage", y="MAF", hue="Biological_group",
                showfliers=False, ax=ax)
    sns.stripplot(data=long, x="Stage", y="MAF", hue="Biological_group",
                  dodge=True, color="black", alpha=.45, size=3, ax=ax, legend=False)
    ax.set(title="Variant allele fraction across tissue stages",
           xlabel="Tissue stage", ylabel="Mutant allele fraction (MAF)")
    fig.tight_layout()
    path = output_dir / "maf_by_stage.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    module = axis.dropna(subset=["Module"])
    if not module.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=module, x="Delta_MAF", y="Module", showfliers=False, ax=ax)
        sns.stripplot(data=module, x="Delta_MAF", y="Module",
                      color="black", alpha=.55, size=3, ax=ax)
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set(title="Change in MAF from actinic keratosis to cSCC",
               xlabel="ΔMAF (cSCC − AK)", ylabel="Pathway module")
        fig.tight_layout()
        path = output_dir / "delta_maf_by_module.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths
