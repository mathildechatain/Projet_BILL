#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import os                                                                                   # importe les packages nécessaires pour le script

argv = sys.argv                                                                             # crée une variable "argv" qui contient les arguments du script
if len(argv) < 2:                                                                           # conditionnelle: si la taille de l'argument est inférieure à 2
    print("Usage: python ReadVCFsnp.py <input_folder_vcf>")                                 # affichage d'un message d'erreur
    sys.exit(1)                                                                             # sortie du script

vcf_folder = argv[1]                                                                        # la variable "vcf_folder" reçoit l'information contenu à la position 1 de argv. Cette information est un chemin de fichier contenant les fichiers *.vcf
List_QUAL = []                                                                              # une liste vide est crée "List_QUAL"

for filename in os.listdir(vcf_folder):                                                     # crée une variable "filename" qui prendra successivement chaque nom des fichiers dans le dossier fourni en argument
    if filename.endswith(".snp.vcf"):                                                       # conditionnelle: si le nom du fichier fini par l'extension ".snp.vcf"
        filepath = os.path.join(vcf_folder, filename)                                       # on récupère le chemin du fichier dans la variable "filepath"
        with open(filepath) as f:                                                           # crée l'objet "f" qui est le fichier ouvert
            for line in f:                                                                  # la variable line parcoura successivement chaque ligne du fichier
                if line.startswith("#"):                                                    # conditionnelle: si la ligne commence par "#"
                    continue                                                                # ignorer la ligne
                fields = line.strip().split("\t")                                           # la variable "fields" est une liste qui contiendra toutes les colonnes d'une ligne
                try:                                                                        #
                    qual = float(fields[5])                                                 # la variable "qual" contiendra la valeur de qualité situé à la 5eme position de "fields"
                    List_QUAL.append(qual)                                                  # La valeur de qual sera ajouté à la liste "List_QUAL"
                except ValueError:                                                          # vérifie qu'il n'y a pas d'erreur
                    continue  # ignorer si vide ou invalide                                 # si oui alors continue

print(f"Nombre de variants : {len(List_QUAL)}")                                             # affiche dans le terminal le nombre de variants
print(f"QUAL: {List_QUAL}")                                                                 # affiche la liste des qualités

# Sauvegarde dans un fichier
script_dir = os.path.dirname(os.path.abspath(__file__))                                     # crée une variable contenant le chemin où ce trouve le script
list_qual_dir = os.path.join(script_dir, "List_QUAL")                                       # création du chemin menant à un dossier "List_QUAL"
os.makedirs(list_qual_dir, exist_ok=True)                                                   # création du dossier "List_QUAL" s'il n'existe pas déjà
output_file_support = os.path.join(list_qual_dir, "List_QUAL_all.txt")                      # création du chemin menant à un fichier texte "List_QUAL_all.txt"
with open(output_file_support, "w") as out:                                                 # ouverture du fichier en mode écriture (création du fichier lors de l'ouverture)
    out.write(",".join(str(q) for q in List_QUAL))                                          # réécrit la liste de qualité dans le fichier texte avec les valeurs séparées par des virgules


if List_QUAL:                                                                               # Si la liste est vide, on ne fait pas le plot
    plt.figure(figsize=(10,6))                                                              # configuration de la taille du plot
    plt.hist(List_QUAL, bins=50, color='skyblue', edgecolor='black')                        # plot de type histogramme prennant en valeurs les qualité par fenêtre de 50 et option couleur
    plt.title("Distribution des valeurs QUAL des polymorphismes nucléotidiques simples")    # donne le titre du graphique
    plt.xlabel("QUAL")                                                                      # nomme l'axe x
    plt.ylabel("Nombre de variants")                                                        # nomme l'axe y
    plt.grid(axis='y', alpha=0.75)                                                          # affiche le quadrillage sur le graphique

# Chemin vers le dossier plot
repo_root = os.path.dirname(script_dir)                                                     # crée une variable contenant le chemin qui mène au dossier en amont de celui contenant le script
plot_folder = os.path.join(repo_root, "Plots", "QUAL_distribution")                         # création du chemin menant à des dossiers imbriqué "Plots" contenant "QUAL_distribution"
os.makedirs(plot_folder, exist_ok=True)                                                     # crée le dossier s'il n'existe pas

plot_path = os.path.join(plot_folder, "Distribution_QUAL_variants_snp.png")                 # création du chemin pour l'image du graphique
plt.savefig(plot_path)                                                                      # sauvegarde du graphique au format png au lieu indiqué
print(f"Plot enregistre dans {plot_path}")                                                  # affiche que le plot a été sauvegaré dans le dossier indiqué
plt.show()                                                                                  # ouvre une fenetre qui affiche le graphique créé
