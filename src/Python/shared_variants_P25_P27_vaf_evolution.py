#!/usr/bin/env python3
from pathlib import Path  
import matplotlib.pyplot as plt  
import pandas as pd  


# Définition des chemins (dossier du script et racine projet)
SCRIPT_DIR = Path(__file__).resolve().parent  
PROJECT_ROOT = SCRIPT_DIR.parent.parent  

# Fichier d'entrée : table des variants SV annotés
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"
# Fichier de sortie csv
OUT_SUMMARY = PROJECT_ROOT / "Files" / "plot_files" / "shared_variants_P25_P27_vaf_evolution_summary.csv"
# Fichier de sortie : graphique dans Plots
OUT_PLOT = PROJECT_ROOT / "Plots" / "shared_variants_P25_P27_vaf_evolution.png"

# Générations étudiées
GEN_ORDER = ["P25", "P27"]
# Colonnes définissant un variant unique
VARIANT_KEYS = ["chrom", "position", "ref", "alt", "svtype", "svlen"]
# Couleurs pour le graphique
BOX_COLOR = "#F6C177"  # couleur des boxplots
POINT_COLOR = "#C05621"  # couleur des points (échantillons)


def main() -> None:
    # Lecture du fichier CSV
    df = pd.read_csv(INPUT_CSV)
    # Conversion des colonnes en numérique
    df["vaf"] = pd.to_numeric(df["vaf"], errors="coerce")  # fréquence allélique
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    # Filtrage des générations d'intérêt
    df = df[df["generation"].isin(GEN_ORDER)].copy()
    # Filtrage des echantillons valides (1 à 10)
    df = df[df["replicate"].between(1, 10)].copy()

    # Identification des variants partagés entre P25 et P27
    shared_variants = (
        df[df["generation"] == "P25"][VARIANT_KEYS]
        .drop_duplicates()  # variants uniques en P25
        .merge(
            df[df["generation"] == "P27"][VARIANT_KEYS].drop_duplicates(),
            on=VARIANT_KEYS,  # intersection avec P27
        )
    )



    # Liste des échantillons observés
    observed_samples = df[["generation", "replicate", "sample"]].drop_duplicates()
    # Produit cartésien : tous les échantillons x tous les variants partagés
    sample_variant_grid = observed_samples.merge(shared_variants, how="cross")
    # Association des VAF observées 
    sample_variant_vaf = sample_variant_grid.merge(
        df[["generation", "replicate", "sample", "vaf"] + VARIANT_KEYS],
        on=["generation", "replicate", "sample"] + VARIANT_KEYS,
        how="left",
    )

    # Remplacement des valeurs manquantes (variants absents) par 0
    sample_variant_vaf["vaf"] = sample_variant_vaf["vaf"].fillna(0.0)

    # Résumé par échantillon
    sample_summary = (
        sample_variant_vaf.groupby(["generation", "replicate", "sample"], as_index=False)
        .agg(
            mean_vaf=("vaf", "mean"),  # moyenne des VAF
            median_vaf=("vaf", "median"),  # médiane
            detected_shared_variants=("vaf", lambda s: int((s > 0).sum())),  # variants détectés
        )
    )

    # Nombre total de variants partagés
    sample_summary["total_shared_variants"] = len(shared_variants)

    # Conversion en catégories ordonnées
    sample_summary["generation"] = pd.Categorical(
        sample_summary["generation"], categories=GEN_ORDER, ordered=True
    )

    # Tri des résultats
    sample_summary = sample_summary.sort_values(["generation", "replicate"])

    # Sauvegarde du résumé
    sample_summary.to_csv(OUT_SUMMARY, index=False)

   
    box_data = [
        sample_summary.loc[sample_summary["generation"] == generation, "mean_vaf"].tolist()
        for generation in GEN_ORDER
    ]
    fig, ax = plt.subplots(figsize=(8.5, 5))

    # Positions des boxplots
    positions = list(range(1, len(GEN_ORDER) + 1))

    # Création des boxplots
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

    # Ajout des points individuels (chaque échantillon)
    for pos, generation in zip(positions, GEN_ORDER):
        gen_values = sample_summary.loc[
            sample_summary["generation"] == generation, ["replicate", "mean_vaf"]
        ]
        if gen_values.empty:
            continue

        #  éviter la superposition des points
        offsets = pd.Series(range(len(gen_values)), index=gen_values.index, dtype=float)
        if len(gen_values) > 1:
            offsets = (offsets - (len(gen_values) - 1) / 2) * 0.05
        else:
            offsets = offsets * 0

        # Scatter plot des valeurs individuelles
        ax.scatter(
            pos + offsets,
            gen_values["mean_vaf"],
            s=35,
            color=POINT_COLOR,
            alpha=0.85,
            zorder=3,
        )

    # Labels du graphique
    ax.set_xticks(positions)
    ax.set_xticklabels(GEN_ORDER)
    ax.set_xlabel("Generation")
    ax.set_ylabel("VAF moyenne des variants partages P25/P27")
    ax.set_title("Evolution de la VAF des variants partages entre P25 et P27")
    ax.set_ylim(bottom=0)

    ax.grid(True, axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()

    # Sauvegarde du graphique
    fig.savefig(OUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Messages affichés
    print("Sauvegardé:")
    print("-", OUT_SUMMARY)
    print("-", OUT_PLOT)
    print("-", f"Shared variants used: {len(shared_variants)}")

if __name__ == "__main__":
    main()