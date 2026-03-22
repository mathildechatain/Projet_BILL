#!/usr/bin/env python3
from pathlib import Path 
import pandas as pd  


# Définition des chemins de dossier du script et racine du projet
SCRIPT_DIR = Path(__file__).resolve().parent  
PROJECT_ROOT = SCRIPT_DIR.parent.parent  

# Fichier d'entrée contenant les variants SV annotés
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"
# Fichier de sortie 1
OUT_COUNTS_LONG = PROJECT_ROOT / "Files" / "plot_files" / "sv_by_orf_generation_summary.csv"
# Fichier de sortie 2
OUT_COUNTS_WIDE = PROJECT_ROOT / "Files" / "plot_files" / "sv_by_orf_generation_wide.csv"
# Ordre des générations
GEN_ORDER = ["P15", "P25", "P27", "P30", "P50", "P65"]


# Fonction pour séparer les ORFs multiples
def split_orfs(value: object) -> list[str]:
    if pd.isna(value) or str(value).strip() == "":
        return []  # aucun ORF
    # séparation par ";" et nettoyage des espaces
    return [part.strip() for part in str(value).split(";") if part.strip()]


def main() -> None:
    # Lecture du fichier CSV
    df = pd.read_csv(INPUT_CSV)
    # Conversion d'echantillon en numérique
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")
    # Filtrage des echantillons valides (1 à 10)
    df = df[df["replicate"].between(1, 10)].copy()
    # Filtrage des générations d'intérêt
    df = df[df["generation"].isin(GEN_ORDER)].copy()
    # Création d'un identifiant unique pour chaque variant
    df["variant_id"] = df[["chrom", "position", "svtype", "ref", "alt", "svlen"]].astype(str).agg("|".join, axis=1)
    # Séparation des ORFs multiples en listes
    df["orf"] = df["orf"].apply(split_orfs)
    # une ligne par ORF
    df = df.explode("orf")
    # Suppression des ORFs vides
    df = df[df["orf"].notna() & (df["orf"] != "")].copy()
    # Suppression des variants intergéniques
    df = df[df["orf"] != "Intergenic"].copy()

    # Comptage des variants uniques par ORF et génération
    counts = (
        df.groupby(["orf", "generation"], as_index=False)["variant_id"]
        .nunique()  # nombre de variants uniques
        .rename(columns={"variant_id": "n_unique_sv"})
    )

    # Conversion de la génération en facteur ordonné
    counts["generation"] = pd.Categorical(
        counts["generation"],
        categories=GEN_ORDER,
        ordered=True
    )

    # Tri des résultats
    counts = counts.sort_values(["orf", "generation"])

    # Sauvegarde au format long
    counts.to_csv(OUT_COUNTS_LONG, index=False)

    
    # Transformation en format large pour visualiser
    wide = (
        counts.pivot(index="orf", columns="generation", values="n_unique_sv")
        .reindex(columns=GEN_ORDER)  # garantir l'ordre des générations
        .fillna(0)  # remplacer les valeurs manquantes par 0
        .astype(int)
    )

    # Calcul du total de SV par ORF
    wide["total_sv"] = wide.sum(axis=1)
    # Tri décroissant selon le total
    wide = wide.sort_values("total_sv", ascending=False)
    # Sauvegarde au format large
    wide.to_csv(OUT_COUNTS_WIDE)

    # Messages de confirmation
    print("Sauvegardé:")
    print("-", str(OUT_COUNTS_LONG))
    print("-", str(OUT_COUNTS_WIDE))

if __name__ == "__main__":
    main()