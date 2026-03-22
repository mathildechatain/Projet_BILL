#!/usr/bin/env python3                                                                                    //pour avoir explication du code voir les commentaires de QUAL_distribution_snp
import sys
import matplotlib.pyplot as plt
import os

argv = sys.argv
if len(argv) < 2:
    print("Usage: python ReadVCFsnp.py <input_folder_vcf>")
    sys.exit(1)

vcf_folder = argv[1]
List_QUAL = []

for filename in os.listdir(vcf_folder):
    if filename.endswith(".snp.vcf") and (filename.startswith("P25") or filename.startswith("P27")):    //ne selectionne que les fichiers des générations P25 et P27
        filepath = os.path.join(vcf_folder, filename)
        with open(filepath) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                try:
                    qual = float(fields[5])
                    List_QUAL.append(qual)
                except ValueError:
                    continue  # ignorer si vide ou invalide

print(f"Nombre de variants : {len(List_QUAL)}")
print(f"QUAL: {List_QUAL}")

# Sauvegarde dans fichier
script_dir = os.path.dirname(os.path.abspath(__file__))
list_qual_dir = os.path.join(script_dir, "List_QUAL")
os.makedirs(list_qual_dir, exist_ok=True)
output_file_support = os.path.join(list_qual_dir, "List_QUAL_P25_P27.txt")
with open(output_file_support, "w") as out:
    out.write(",".join(str(q) for q in List_QUAL))

# Si la liste est vide, on ne fait pas le plot
if List_QUAL:
    plt.figure(figsize=(10,6))
    plt.hist(List_QUAL, bins=50, color='skyblue', edgecolor='black')
    plt.title("Distribution des valeurs QUAL des polymorphismes nucléotidiques simples pour P25 et P27")
    plt.xlabel("QUAL")
    plt.ylabel("Nombre de variants")
    plt.grid(axis='y', alpha=0.75)

# Chemin vers le dossier plot
repo_root = os.path.dirname(script_dir)
plot_folder = os.path.join(repo_root, "Plots", "QUAL_distribution")
os.makedirs(plot_folder, exist_ok=True)

plot_path = os.path.join(plot_folder, "Distribution_QUAL_variants_snp_P25_P27.png")
plt.savefig(plot_path)
print(f"Plot enregistre dans {plot_path}")
plt.show()
