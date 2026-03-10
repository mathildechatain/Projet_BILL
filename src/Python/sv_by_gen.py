#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt


INPUT_CSV = "../../Files/plot_files/summary_table_sv_with_ORF.csv"
OUT_COUNTS = "../../Files/plot_files/sv_per_generation_summary.csv"
OUT_PLOT = "../../Files/plot_files/sv_per_generation.png"

GEN_ORDER = ["P15", "P25", "P27", "P30", "P50", "P65"]


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    df = df[df["replicate"].between(1, 10)].copy()
    df = df[df["generation"].isin(GEN_ORDER)].copy()

    # Number of SV per sample (replicate) and generation
    sample_counts = (
        df.groupby(["generation", "replicate"], as_index=False)
        .size()
        .rename(columns={"size": "n_sv"})
    )

    # Mean evolution per generation
    summary = (
        sample_counts.groupby("generation", as_index=False)["n_sv"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "mean_n_sv", "std": "sd_n_sv", "count": "n_samples"})
    )
    summary["sem_n_sv"] = summary["sd_n_sv"] / (summary["n_samples"] ** 0.5)
    summary["generation"] = pd.Categorical(summary["generation"], categories=GEN_ORDER, ordered=True)
    summary = summary.sort_values("generation")
    summary.to_csv(OUT_COUNTS, index=False)

    plt.figure(figsize=(8, 5))
    x = list(range(len(summary)))
    plt.plot(x, summary["mean_n_sv"], marker="o", linewidth=2, color="#2B6CB0")
    plt.errorbar(
        x,
        summary["mean_n_sv"],
        yerr=summary["sem_n_sv"],
        fmt="none",
        ecolor="#2B6CB0",
        elinewidth=1.2,
        capsize=4,
    )
    plt.xticks(x, summary["generation"])
    plt.xlabel("Generation")
    plt.ylabel("Nombre moyen de SV par echantillon")
    plt.title("Evolution du nombre moyen de mutations SV par generation")
    plt.grid(True, axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:")
    print("-", OUT_COUNTS)
    print("-", OUT_PLOT)


if __name__ == "__main__":
    main()
