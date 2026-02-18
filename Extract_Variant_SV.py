import os
import glob
from pathlib import Path
import sys

# Vérifie qu'un argument est fourni
if len(sys.argv) < 2:
    print("Usage: python filter_snp_vcf.py <dossier_vcf>") #sinon affiche messaege d'erreur
    sys.exit(1)

vcf_folder = sys.argv[1]  # Le dossier passé en argument
output_folder = os.path.join(vcf_folder, "vcf_filtered")
Path(output_folder).mkdir(exist_ok=True)

# Fonction pour filtrer les SV sur notre critère QUAL > 10, qui prend en entrée un fichier VCF et un fichier de sortie
def filter_SV_vcf(input_file, output_file):

    # Ouvrir le fichier d'entrée et le fichier de sortie
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Parcourir chaque ligne du fichier d'entrée
        for line in infile:
            if line.startswith("#"):
                # Garder toutes les lignes d'en-tête
                outfile.write(line)
            else:
                # Split la ligne en colonnes pour accéder à la colonne QUAL (index 5)
                cols = line.strip().split("\t")
                try:
                    # Accéder à la colonne QUAL (index 5) et vérifier si elle est supérieure à 10
                    DV = float(cols[9].split(":")[3])  # DV est la 4ème valeur du champ FORMAT
                    if DV >= 200:
                        # Si DV est supérieur à 10, écrire la ligne dans le fichier de sortie
                        outfile.write(line)
                except ValueError:
                    # Ignore si DV n'est pas un nombre
                    continue

# Traitement de tous les fichiers VCF
# Utiliser glob pour trouver tous les fichiers VCF dans le dossier spécifié
vcf_files = glob.glob(os.path.join(vcf_folder, "*sniffles.vcf"))
# Parcourir chaque fichier VCF trouvé et appliquer la fonction de filtrage
for vcf_file in vcf_files:
    # Générer le nom du fichier de sortie en ajoutant "_SV_filtered" au nom de base du fichier d'entrée
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    # Créer le chemin complet du fichier de sortie
    output_file = os.path.join(output_folder, f"{base_name}_SV_filtered.vcf")
    # Appeler la fonction de filtrage pour chaque fichier VCF
    filter_SV_vcf(vcf_file, output_file)

print("Filtrage SV terminé !")

