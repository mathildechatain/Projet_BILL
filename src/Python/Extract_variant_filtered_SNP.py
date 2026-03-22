import os
import glob
from pathlib import Path
import sys                                                                          # importe les packages nécessaires pour le script


if len(sys.argv) < 2:                                                               # conditionnelle: si la taille de l'argument est inférieure à 2
    print("Usage: python filter_snp_vcf.py <dossier_vcf>")                          # affichage d'un message d'erreur
    sys.exit(1)                                                                     # sortie du script


vcf_folder = sys.argv[1]                                                            # la variable "vcf_folder" reçoit l'information contenu à la position 1 de argv. Cette information est un chemin de fichier contenant les fichiers *.vcf
output_folder = os.path.join("../../Files", "vcf_filtered_snp")                     # crée une variable contenant le chemin où ce trouve le dossier "vcf_filtered_snp"
Path(output_folder).mkdir(exist_ok=True)                                            # création du dossier, s'il n'existe pas déjà

def filter_snp_vcf(input_file, output_file):                                        # Fonction pour filtrer les SNP en fonction de notre critère QUAL > 20, qui prend en entrée un fichier VCF et un fichier de sortie

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:        # ouvre le fichier VCF d'entrée en lecture et le fichier de sortie en écriture
        for line in infile:                                                         # la variable "line" parcoura successivement chaque ligne du fichier d'entrée
            if line.startswith("#"):                                                # conditionnelle si la ligne commence par "#"
                outfile.write(line)                                                 # Conserve toutes les lignes d'en-tête (noms de colonnes) inchangées dans le fichier de sortie
            else:
                cols = line.strip().split("\t")                                     # la variable "fields" est une liste qui contiendra toutes les colonnes d'une ligne en utilisant la tabulation comme délimiteur
                try:
                    qual = float(cols[5])                                           # la variable "qual" contiendra la valeur de qualité situé à la 5eme position de "fields"
                    if qual > 20:                                                   # conditionnelle: si la qualité est supérieure à 20
                        outfile.write(line)                                         # écrit la ligne dans le fichier de sortie
                except ValueError:                                                  # vérifie qu'il n'y a pas d'erreur
                    continue                                                        # si oui alors continue



vcf_files = glob.glob(os.path.join(vcf_folder, "*snp.vcf"))                         # Utilise glob pour trouver tous les fichiers VCF dans le dossier spécifié qui se terminent par "snp.vcf", ce qui correspond aux variants SNP
for vcf_file in vcf_files:                                                          # Itère sur chaque fichier VCF trouvé et applique la fonction de filtrage
    base_name = os.path.basename(vcf_file).replace(".vcf", "")                      # Génère le nom du fichier de sortie en ajoutant "_SNP_filtered" au nom de base du fichier d'entrée
    output_file = os.path.join(output_folder, f"{base_name}_SNP_filtered.vcf")      # Crée le chemin complet pour le fichier de sortie
    filter_snp_vcf(vcf_file, output_file)                                           # Appelle la fonction de filtrage pour chaque fichier VCF


print("SNP filtering completed!")                                                   # Affiche un message pour indiquer que le processus de filtrage est terminé

