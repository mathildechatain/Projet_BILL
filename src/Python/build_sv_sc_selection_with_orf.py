#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"
OUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "sv_sc_turnover.csv"

VARIANT_KEYS = ["chrom", "position", "ref", "alt", "svtype", "svlen"]
PHASES = [
    ("control_P15_P25", "P15", "P25"),
    ("stressed_P25_P27", "P25", "P27"),
    ("control_P30_P50", "P30", "P50"),
]


def lineage_group(replicate: int) -> str:
    return "lineages_1_5" if replicate <= 5 else "lineages_6_10"


def normalize_orfs(series: pd.Series) -> str:
    values = set()
    for entry in series.dropna():
        for item in str(entry).split(";"):
            item = item.strip()
            if item:
                values.add(item)
    return ";".join(sorted(values))


def variants_for_sample(df: pd.DataFrame, generation: str, replicate: int) -> pd.DataFrame:
    subset = df[(df["generation"] == generation) & (df["replicate"] == replicate)].copy()
    subset["variant_id"] = pd.Series(dtype=str)
    if subset.empty:
        return subset
    subset["variant_id"] = subset[VARIANT_KEYS].astype(str).agg("|".join, axis=1)
    return subset.drop_duplicates("variant_id")


def category_row(phase: str, g_start: str, g_end: str, replicate: int, category: str, variants: pd.DataFrame) -> dict:
    return {
        "phase": phase,
        "g_start": g_start,
        "g_end": g_end,
        "replicate": replicate,
        "lineage_group": lineage_group(replicate),
        "category": category,
        "n_variants": len(variants),
        "orfs": normalize_orfs(variants["orf"]) if not variants.empty else "",
    }


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce").astype("Int64")

    rows = []
    for phase, g_start, g_end in PHASES:
        for replicate in range(1, 11):
            start_df = variants_for_sample(df, g_start, replicate)
            end_df = variants_for_sample(df, g_end, replicate)

            start_ids = set(start_df["variant_id"]) if not start_df.empty else set()
            end_ids = set(end_df["variant_id"]) if not end_df.empty else set()

            shared_ids = start_ids & end_ids
            gained_ids = end_ids - start_ids
            lost_ids = start_ids - end_ids

            rows.append(
                category_row(
                    phase,
                    g_start,
                    g_end,
                    replicate,
                    "shared",
                    end_df[end_df["variant_id"].isin(shared_ids)],
                )
            )
            rows.append(
                category_row(
                    phase,
                    g_start,
                    g_end,
                    replicate,
                    "gained_end",
                    end_df[end_df["variant_id"].isin(gained_ids)],
                )
            )
            rows.append(
                category_row(
                    phase,
                    g_start,
                    g_end,
                    replicate,
                    "lost_from_start",
                    start_df[start_df["variant_id"].isin(lost_ids)],
                )
            )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_CSV, index=False)

    print("Saved:")
    print("-", OUT_CSV)


if __name__ == "__main__":
    main()
