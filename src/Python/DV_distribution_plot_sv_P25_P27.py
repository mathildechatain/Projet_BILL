#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import os

#Objectif du script :
#Créer une liste de toutes les valeurs DV des variants structuraux pour P25 et P27 à partir des fichiers VCF sniffles de type SV
#puis faire un histogramme de la distribution de ces valeurs DV

argv=sys.argv
if len(argv)<2: #si l'argument a une taille supérieure à 2
    print("Usage: python ReadVCFsv.py <input.vcf>")
    sys.exit(1) 

vcf_folder=argv[1]
List_DV=[]

for filename in os.listdir(vcf_folder):
         if filename.endswith(".vcf") and (filename.startswith("P27") or filename.startswith("P25")) and "sniffles" in filename:
            filepath = os.path.join(vcf_folder, filename)
            with open(filepath) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields=line.strip().split("\t")
                    chrom=fields[0]
                    pos=fields[1]
                    id=fields[2]
                    filter=fields[6]
                    format=fields[9].split(":")
                    DV=format[3]  # DV est la 4ème valeur du champ FORMAT
                    List_DV.append(DV)

print(f"DV : {List_DV}")

script_dir = os.path.dirname(os.path.abspath(__file__))
list_dv_dir = os.path.join(script_dir, "List_DV")
os.makedirs(list_dv_dir, exist_ok=True)
output_file_support = os.path.join(list_dv_dir, "List_DV_P25_P27.txt")
with open(output_file_support, "w") as out:
    out.write(",".join(str(q) for q in List_DV))

print(f"Liste DV enregistrée dans {output_file_support} (valeurs séparées par des virgules)")
if not List_DV:
    print("Aucune valeur DV trouvee pour P25/P27 avec ce dossier d'entree.")
    sys.exit(0)
with open(output_file_support, "r") as f:
    line = f.readline()  # tout est sur une seule ligne
    DV_values = [float(x) for x in line.strip().split(",")]  # convertir en float

plt.figure(figsize=(10,6))
plt.hist(DV_values, bins=50, color='skyblue', edgecolor='black')
plt.title("Distribution des valeurs DV des variants structuraux pour P25 et P27")
plt.xlabel("DV ( Depth Variant)")
plt.ylabel("Nombre de variants")
plt.grid(axis='y', alpha=0.75)

# Chemin vers la racine du projet
repo_root = os.path.dirname(script_dir)

# Chemin vers le dossier Plots/DV_distribution
plot_folder = os.path.join(repo_root, "Plots", "DV_distribution")

# 🔹 Créer le dossier s'il n'existe pas
os.makedirs(plot_folder, exist_ok=True)

# Chemin complet du fichier
plot_path = os.path.join(plot_folder, "Distribution_DV_variants_SV_P25_P27.png")

# Enregistrer le plot
plt.savefig(plot_path)
print(f"Plot enregistré dans {plot_path}")
plt.show()

