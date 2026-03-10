#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"
OUT_SUMMARY = PROJECT_ROOT / "Files" / "plot_files" / "shared_variants_P25_P27_vaf_evolution_summary.csv"
OUT_PLOT = PROJECT_ROOT / "Plots" / "shared_variants_P25_P27_vaf_evolution.png"

GEN_ORDER = ["P25", "P27"]
VARIANT_KEYS = ["chrom", "position", "ref", "alt", "svtype", "svlen"]
BOX_COLOR = "#F6C177"
POINT_COLOR = "#C05621"


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    df["vaf"] = pd.to_numeric(df["vaf"], errors="coerce")
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    df = df[df["generation"].isin(GEN_ORDER)].copy()
    df = df[df["replicate"].between(1, 10)].copy()

    shared_variants = (
        df[df["generation"] == "P25"][VARIANT_KEYS]
        .drop_duplicates()
        .merge(
            df[df["generation"] == "P27"][VARIANT_KEYS].drop_duplicates(),
            on=VARIANT_KEYS,
        )
    )

    observed_samples = df[["generation", "replicate", "sample"]].drop_duplicates()
    sample_variant_grid = observed_samples.merge(shared_variants, how="cross")
    sample_variant_vaf = sample_variant_grid.merge(
        df[["generation", "replicate", "sample", "vaf"] + VARIANT_KEYS],
        on=["generation", "replicate", "sample"] + VARIANT_KEYS,
        how="left",
    )
    sample_variant_vaf["vaf"] = sample_variant_vaf["vaf"].fillna(0.0)

    sample_summary = (
        sample_variant_vaf.groupby(["generation", "replicate", "sample"], as_index=False)
        .agg(
            mean_vaf=("vaf", "mean"),
            median_vaf=("vaf", "median"),
            detected_shared_variants=("vaf", lambda s: int((s > 0).sum())),
        )
    )
    sample_summary["total_shared_variants"] = len(shared_variants)
    sample_summary["generation"] = pd.Categorical(
        sample_summary["generation"], categories=GEN_ORDER, ordered=True
    )
    sample_summary = sample_summary.sort_values(["generation", "replicate"])
    sample_summary.to_csv(OUT_SUMMARY, index=False)

    box_data = [
        sample_summary.loc[sample_summary["generation"] == generation, "mean_vaf"].tolist()
        for generation in GEN_ORDER
    ]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    positions = list(range(1, len(GEN_ORDER) + 1))
    ax.boxplot(
        box_data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": BOX_COLOR, "edgecolor": POINT_COLOR, "linewidth": 1.2},
        medianprops={"color": "#1A202C", "linewidth": 1.5},
        whiskerprops={"color": POINT_COLOR, "linewidth": 1.2},
        capprops={"color": POINT_COLOR, "linewidth": 1.2},
    )

    for pos, generation in zip(positions, GEN_ORDER):
        gen_values = sample_summary.loc[
            sample_summary["generation"] == generation, ["replicate", "mean_vaf"]
        ]
        if gen_values.empty:
            continue
        offsets = pd.Series(range(len(gen_values)), index=gen_values.index, dtype=float)
        if len(gen_values) > 1:
            offsets = (offsets - (len(gen_values) - 1) / 2) * 0.05
        else:
            offsets = offsets * 0
        ax.scatter(
            pos + offsets,
            gen_values["mean_vaf"],
            s=35,
            color=POINT_COLOR,
            alpha=0.85,
            zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(GEN_ORDER)
    ax.set_xlabel("Generation")
    ax.set_ylabel("VAF moyenne des variants partages P25/P27")
    ax.set_title("Evolution de la VAF des variants partages entre P25 et P27")
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved:")
    print("-", OUT_SUMMARY)
    print("-", OUT_PLOT)
    print("-", f"Shared variants used: {len(shared_variants)}")


if __name__ == "__main__":
    main()
