import os
import glob
from pathlib import Path
import sys

# vérifie que l'utilisateur a fourni le dossier contenant les fichiers VCF en argument
if len(sys.argv) < 2:
    print("Usage: python filter_snp_vcf.py <dossier_vcf>") # si ce n'est pas le cas, affiche un message d'utilisation et quitte le programme avec un code d'erreur    
    sys.exit(1)

# le dossier contenant les fichiers VCF est passé en argument
vcf_folder = sys.argv[1]
# Crée le dossier de sortie, nommé vcf_filtered, s'il n'existe pas
output_folder = os.path.join("../Files", "vcf_filtered_snp")
Path(output_folder).mkdir(exist_ok=True)

# Fonction pour filtrer les SNP en fonction de notre critère QUAL > 10, qui prend en entrée un fichier VCF et un fichier de sortie
def filter_snp_vcf(input_file, output_file):

    # ouvre le fichier VCF d'entrée en lecture et le fichier de sortie en écriture
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Itère sur chaque ligne du fichier d'entrée
        for line in infile:
            if line.startswith("#"):
                # Conserve toutes les lignes d'en-tête (noms de colonnes) inchangées dans le fichier de sortie
                outfile.write(line)
            else:
                # Divise la ligne en colonnes en utilisant la tabulation comme délimiteur
                cols = line.strip().split("\t")
                try:
                    # Accède à la colonne QUAL (index 5) et vérifie si elle est supérieure à 10
                    qual = float(cols[5])
                    if qual > 20:
                        # Si QUAL est supérieure à 10, écrit la ligne dans le fichier de sortie
                        outfile.write(line)
                except ValueError:
                    continue

# Traitement de tous les fichiers VCF
# Utilise glob pour trouver tous les fichiers VCF dans le dossier spécifié qui se terminent par "snp.vcf", ce qui correspond aux variants SNP
vcf_files = glob.glob(os.path.join(vcf_folder, "*snp.vcf"))
# Itère sur chaque fichier VCF trouvé et applique la fonction de filtrage
for vcf_file in vcf_files:
    # Génère le nom du fichier de sortie en ajoutant "_SNP_filtered" au nom de base du fichier d'entrée
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    # Crée le chemin complet pour le fichier de sortie
    output_file = os.path.join(output_folder, f"{base_name}_SNP_filtered.vcf")
    # Appelle la fonction de filtrage pour chaque fichier VCF
    filter_snp_vcf(vcf_file, output_file)

# Affiche un message pour indiquer que le processus de filtrage est terminé
print("SNP filtering completed!")

