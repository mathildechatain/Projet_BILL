#!/usr/bin/env python3
from pathlib import Path 
import pandas as pd  

# Définition des chemins principaux ( dossier script et racine projet )
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent 

# Fichier d'entrée contenant les variants SV annotés
INPUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "summary_table_sv_with_ORF.csv"
# Fichier de sortie CSV utile pour plot R
OUT_CSV = PROJECT_ROOT / "Files" / "plot_files" / "sv_sc_turnover.csv"

# Colonnes utilisées pour définir un identifiant unique de variant
VARIANT_KEYS = ["chrom", "position", "ref", "alt", "svtype", "svlen"]
# Définition des passages (comparaisons entre générations)
PASSAGES = [
    ("control_P15_P25", "P15", "P25"),
    ("stressed_P25_P27", "P25", "P27"),
    ("control_P30_P50", "P30", "P50"),
]


# Attribution du groupe d'echantillons selon le replicate
def lineage_group(replicate: int) -> str:
    return "lineages_1_5" if replicate <= 5 else "lineages_6_10"


# Normalisation des ORFs : supprime les doublons et trie les valeurs
def normalize_orfs(series: pd.Series) -> str:
    values = set()  # Utilisation d'un set pour éviter les doublons
    for entry in series.dropna():  
        for item in str(entry).split(";"):  # Sépare les ORFs multiples
            item = item.strip()
            if item:
                values.add(item)  
    return ";".join(sorted(values)) #chaîne triée


# Extraction des variants pour une génération et un échantillon donnés
def variants_for_sample(df: pd.DataFrame, generation: str, replicate: int) -> pd.DataFrame:
    # Filtrage du DataFrame
    subset = df[(df["generation"] == generation) & (df["replicate"] == replicate)].copy()
    # Initialisation de la colonne variant_id
    subset["variant_id"] = pd.Series(dtype=str)
    # Si aucun variant, retourner directement
    if subset.empty:
        return subset

    # Création d'un identifiant unique pour chaque variant
    subset["variant_id"] = subset[VARIANT_KEYS].astype(str).agg("|".join, axis=1)

    # Suppression des doublons
    return subset.drop_duplicates("variant_id")


# Création d'une ligne de résultat pour une catégorie donnée
def category_row(passage: str, g_start: str, g_end: str, replicate: int, category: str, variants: pd.DataFrame) -> dict:
    return {
        "passage": passage,  # Nom du passage
        "g_start": g_start,  # Génération de départ
        "g_end": g_end,  # Génération de fin
        "replicate": replicate,  # Numéro du replicate
        "lineage_group": lineage_group(replicate),  # Groupe de lignées
        "category": category,  # Type de variants (shared / gained / lost)
        "n_variants": len(variants),  # Nombre de variants dans cette catégorie
        # ORFs associés
        "orfs": normalize_orfs(variants["orf"]) if not variants.empty else "",
    }


def main() -> None:
    # Lecture du fichier CSV
    df = pd.read_csv(INPUT_CSV)
    # Conversion de la colonne replicate en entier 
    df["replicate"] = pd.to_numeric(df["replicate"], errors="coerce").astype("Int64")
    rows = []  # Liste pour stocker les résultats

    # Boucle sur chaque passage (comparaison de générations)
    for passage, g_start, g_end in PASSAGES:
        # Boucle sur les replicates de 1 à 10
        for replicate in range(1, 11):
            # Extraction des variants au début et à la fin du passage
            start_df = variants_for_sample(df, g_start, replicate)
            end_df = variants_for_sample(df, g_end, replicate)

            # Création des ensembles d'identifiants de variants
            start_ids = set(start_df["variant_id"]) if not start_df.empty else set()
            end_ids = set(end_df["variant_id"]) if not end_df.empty else set()

            # Comparaison des ensembles
            shared_ids = start_ids & end_ids  # Variants présents aux deux temps
            gained_ids = end_ids - start_ids  # Variants apparus
            lost_ids = start_ids - end_ids  # Variants disparus

            # Ajout des variants partagés
            rows.append(
                category_row(
                    passage,
                    g_start,
                    g_end,
                    replicate,
                    "shared",
                    end_df[end_df["variant_id"].isin(shared_ids)],
                )
            )

            # Ajout des variants gagnés (présents à la fin uniquement)
            rows.append(
                category_row(
                    passage,
                    g_start,
                    g_end,
                    replicate,
                    "gained_end",
                    end_df[end_df["variant_id"].isin(gained_ids)],
                )
            )

            # Ajout des variants perdus (présents au début uniquement)
            rows.append(
                category_row(
                    passage,
                    g_start,
                    g_end,
                    replicate,
                    "lost_from_start",
                    start_df[start_df["variant_id"].isin(lost_ids)],
                )
            )

    # Création du DataFrame 
    out_df = pd.DataFrame(rows)
    # Sauvegarde en CSV
    out_df.to_csv(OUT_CSV, index=False)

    # Message de confirmation
    print("Sauvegardé:")
    print("-", OUT_CSV)

if __name__ == "__main__":
    main()