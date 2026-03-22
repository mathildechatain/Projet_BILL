#!/usr/bin/env python3
from pathlib import Path  # Gestion des chemins de fichiers

import pandas as pd  # Manipulation de données tabulaires


# Définition des chemins principaux ( doossier du script, racine du projet )
SCRIPT_DIR = Path(__file__).resolve().parent 
PROJECT_ROOT = SCRIPT_DIR.parent.parent  

# Fichier CSV d'entrée (table récapitulative des SVs)
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"

# Fichier CSV de sortie dans Files_plot_files
OUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "sv_sc_counts.csv"

# Définition des passages avec les générations associées
PASSAGES = [
    ("control_P15_P25", ["P15", "P25"]),  # passage contrôle avant choc thermique
    ("stressed_P25_P27", ["P25", "P27"]),  # passage de stress P25-P27 (choc thermique en P26)
    ("control_P30_P50", ["P30", "P50"]),  # passage contrôle après choc
]

# Fonction pour déterminer le groupe échantillonnal en fonction du numéro ( 1 à 5 ou 5 à 10 car traitement thermique différent )
def lineage_group(replicate: float) -> str:
    return "lineages_1_5" if replicate <= 5 else "lineages_6_10"


def main() -> None:
    # Lecture du fichier CSV d'entrée
    df = pd.read_csv(INPUT_CSV)

    # Conversion de la colonne "replicate" en numérique ( si valeurs invalides -> NaN)
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce")

    # Compte le nombre de SV (Structural Variants) par génération et pat échantillon
    sample_counts = (
        df.groupby(["generation", "replicate"], as_index=False)
        .size()  
        .rename(columns={"size": "n_SV"})  # renomme la colonne 
    )

    rows = []

    # Boucle sur chaque passage
    for passage, generations in PASSAGES:
        # Filtrage des données pour les générations du passage
        passage_counts = sample_counts[sample_counts["generation"].isin(generations)].copy()

        # Ajout de la colonne passage
        passage_counts["passage"] = passage

        # Attribution du groupe de lignées
        passage_counts["lineage_group"] = passage_counts["replicate"].apply(lineage_group)

        # Sélection des colonnes utiles et stockage temporaire
        rows.append(passage_counts[["passage", "generation", "replicate", "lineage_group", "n_SV"]])

    # Concaténation de tous les passages en un seul DataFrame
    out_df = pd.concat(rows, ignore_index=True)

    # Ordre pour les générations
    generation_order = {"P15": 1, "P25": 2, "P27": 3, "P30": 4, "P50": 5}

    # Ordre pour les passages
    passage_order = {"control_P15_P25": 1, "stressed_P25_P27": 2, "control_P30_P50": 3}

    # Tri des données selon passage, génération et replicate
    out_df = out_df.sort_values(
        by=["passage", "generation", "replicate"],
        key=lambda col: col.map(passage_order if col.name == "passage" else generation_order)
        if col.name in {"passage", "generation"} else col,
    )

    # Sauvegarde de ce dataframe en CSV pour utilisation R
    out_df.to_csv(OUT_CSV, index=False)

    # Message affiché de sauvegarde
    print("Sauvegardé:")
    print("-", OUT_CSV)



if __name__ == "__main__":
    main()
