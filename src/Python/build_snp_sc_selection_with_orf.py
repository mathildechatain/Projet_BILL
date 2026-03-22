#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_snp_with_ORF.csv"
OUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "snp_turnover.csv"

# SNP definition
VARIANT_KEYS = ["chrom", "position", "ref", "alt"]

PHASES = [
    ("control_P15_P25", "P15", "P25"),
    ("stressed_P25_P27", "P25", "P27"),
    ("control_P30_P50", "P30", "P50"),
]


def normalize_orfs(series: pd.Series) -> str:
    """Merge ORF annotations into unique sorted string"""
    values = set()
    for entry in series.dropna():
        for item in str(entry).split(";"):
            item = item.strip()
            if item:
                values.add(item)
    return ";".join(sorted(values))


def variants_for_sample(df: pd.DataFrame, generation: str, replicate: int) -> pd.DataFrame:
    """Extract unique variants for one generation + replicate"""
    subset = df[
        (df["generation"] == generation) &
        (df["replicate"] == replicate)
    ].copy()

    if subset.empty:
        return subset

    subset["variant_id"] = subset[VARIANT_KEYS].astype(str).agg("|".join, axis=1)

    return subset.drop_duplicates("variant_id")


def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    # sécurise replicate
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce").astype("Int64")

    all_rows = []

    for phase, g_start, g_end in PHASES:
        for replicate in range(1, 11):

            start_df = variants_for_sample(df, g_start, replicate)
            end_df = variants_for_sample(df, g_end, replicate)

            start_ids = set(start_df["variant_id"])
            end_ids = set(end_df["variant_id"])

            shared_ids = start_ids & end_ids
            gained_ids = end_ids - start_ids
            lost_ids = start_ids - end_ids

            # Shared
            shared = end_df[end_df["variant_id"].isin(shared_ids)].copy()
            shared["category"] = "shared"
            shared["phase"] = phase
            shared["replicate"] = replicate
            all_rows.append(shared)

            # Gained
            gained = end_df[end_df["variant_id"].isin(gained_ids)].copy()
            gained["category"] = "gained"
            gained["phase"] = phase
            gained["replicate"] = replicate
            all_rows.append(gained)

            # Lost
            lost = start_df[start_df["variant_id"].isin(lost_ids)].copy()
            lost["category"] = "lost"
            lost["phase"] = phase
            lost["replicate"] = replicate
            all_rows.append(lost)

    # concat tous les variants
    result = pd.concat(all_rows, ignore_index=True)

    # agrégation finale
    result = (
        result
        .groupby(["phase", "replicate"] + VARIANT_KEYS + ["category"])
        .agg(
            n_variants=("variant_id", "size"),
            orf=("orf", normalize_orfs)
        )
        .reset_index()
    )

    # sauvegarde
    result.to_csv(OUT_CSV, index=False)

    print("Saved:")
    print("-", OUT_CSV)


if __name__ == "__main__":
    main()
