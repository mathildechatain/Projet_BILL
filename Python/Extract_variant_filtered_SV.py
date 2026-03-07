import os
import glob
import re
from pathlib import Path
import sys

# ouvre le script en vérifiant que l'utilisateur a fourni le dossier contenant les fichiers VCF en argument
if len(sys.argv) < 2:
    print("Usage: python filter_snp_vcf.py <dossier_vcf>") # else print an error message and exit
    sys.exit(1)

# le dossier contenant les fichiers VCF est passé en argument
vcf_folder = sys.argv[1]

# Résout le dossier de sortie à partir de la racine du projet, indépendamment du répertoire d'exécution.
output_folder = os.path.join("..", "Files", "vcf_filtered_sv")
Path(output_folder).mkdir(parents=True, exist_ok=True)

def get_dv_threshold(filename):
    """Retourne le seuil DV spécifique à la génération en fonction du nom de fichier."""
    if filename.startswith("P15"):
        return 200
    if filename.startswith(("P25", "P27")):
        return 200
    if filename.startswith("P30"):
        return 200
    if filename.startswith("P50"):
        return 200
    if filename.startswith("P65"):
        return 200
    return None


# Fonction pour filtrer les variants structuraux en fonction de notre critère DV >= seuil spécifique à la génération, qui prend en entrée un fichier VCF et un fichier de sortie
def filter_SV_vcf(input_file, output_file, dv_threshold):

    # ouvre le fichier VCF d'entrée en lecture et le fichier de sortie en écriture
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # lit le fichier d'entrée ligne par ligne
        for line in infile:
            if line.startswith("#"):
                # conserve toutes les lignes d'en-tête inchangées dans le fichier de sortie (noms de colonnes)
                outfile.write(line)
            else:
                # Divise la ligne en colonnes en utilisant la tabulation comme délimiteur
                cols = line.strip().split("\t")
                try:
                    # Accède à la colonne DV (index 9, 4ème valeur de FORMAT)
                    DV = float(cols[9].split(":")[3])
                    # Filtre la ligne en fonction du seuil DV spécifique à la génération
                    if DV >= dv_threshold:
                        # Si la condition est remplie, écrit la ligne dans le fichier de sortie
                        outfile.write(line)
                except ValueError:
                    continue

# ne garde que les fichiers nommés exactement comme:
# P15-1.trimed1000.sv_sniffles.vcf (numéro génération/réplicat variable)
name_pattern = re.compile(r"^P\d+-\d+\.trimed1000\.sv_sniffles\.vcf$")
vcf_files = [
    p for p in glob.glob(os.path.join(vcf_folder, "*.vcf"))
    if name_pattern.match(os.path.basename(p))
]

# Itère sur chaque fichier VCF et applique la fonction de filtrage
for vcf_file in vcf_files:
    # Génère le nom du fichier de sortie en ajoutant "_SV_filtered" au nom de base du fichier d'entrée
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    dv_threshold = get_dv_threshold(base_name)
    if dv_threshold is None:
        print(f"Skipped (no DV rule): {base_name}")
        continue
    # Crée le chemin complet pour le fichier de sortie
    output_file = os.path.join(output_folder, f"{base_name}_SV_filtered.vcf")
    # Appelle la fonction de filtrage pour chaque fichier VCF
    filter_SV_vcf(vcf_file, output_file, dv_threshold)
    print(f"Filtered {base_name} with DV >= {dv_threshold}")

# Affiche un message pour indiquer que le processus de filtrage est terminé
print("SV filtering completed!")
