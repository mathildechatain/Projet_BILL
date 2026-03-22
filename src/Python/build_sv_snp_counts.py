#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_snp_with_ORF.csv"
OUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "snp_sc_counts.csv"

PHASES = [
    ("control_P15_P25", ["P15", "P25"]),
    ("stressed_P25_P27", ["P25", "P27"]),
    ("control_P30_P50", ["P30", "P50"]),
]


def lineage_group(replicate: float) -> str:
    return "lineages_1_5" if replicate <= 5 else "lineages_6_10"


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")

    sample_counts = (
        df.groupby(["generation", "replicate"], as_index=False)
        .size()
        .rename(columns={"size": "n_SNP"})
    )

    rows = []
    for phase, generations in PHASES:
        phase_counts = sample_counts[sample_counts["generation"].isin(generations)].copy()
        phase_counts["phase"] = phase
        phase_counts["lineage_group"] = phase_counts["replicate"].apply(lineage_group)
        rows.append(phase_counts[["phase", "generation", "replicate", "lineage_group", "n_SNP"]])

    out_df = pd.concat(rows, ignore_index=True)
    generation_order = {"P15": 1, "P25": 2, "P27": 3, "P30": 4, "P50": 5}
    phase_order = {"control_P15_P25": 1, "stressed_P25_P27": 2, "control_P30_P50": 3}
    out_df = out_df.sort_values(
        by=["phase", "generation", "replicate"],
        key=lambda col: col.map(phase_order if col.name == "phase" else generation_order) if col.name in {"phase", "generation"} else col,
    )
    out_df.to_csv(OUT_CSV, index=False)

    print("Saved:")
    print("-", OUT_CSV)


if __name__ == "__main__":
    main()
