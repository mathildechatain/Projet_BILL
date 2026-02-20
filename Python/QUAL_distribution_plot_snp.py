#!/usr/bin/env python3
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
    if filename.endswith(".snp.vcf"):  # assure-toi que c'est exactement ce suffixe
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
output_file = os.path.join(vcf_folder, "QUAL_list.txt")
with open(output_file, "w") as out:
    out.write(",".join(str(q) for q in List_QUAL))
print(f"Liste QUAL enregistrée dans {output_file}")

# Si la liste est vide, on ne fait pas le plot
if List_QUAL:
    plt.figure(figsize=(10,6))
    plt.hist(List_QUAL, bins=50, color='skyblue', edgecolor='black')
    plt.title("Distribution des valeurs QUAL des variants")
    plt.xlabel("QUAL")
    plt.ylabel("Nombre de variants")
    plt.grid(axis='y', alpha=0.75)

    # Chemin vers le dossier plot
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_folder = os.path.join(repo_root, "Plots")
    os.makedirs(plot_folder, exist_ok=True)
    plot_path = os.path.join(plot_folder, "Distribution_QUAL_variants_snp.png")

    plt.savefig(plot_path)
    print(f"Plot enregistré dans {plot_path}")
    plt.show()
else:
    print("Aucun variant QUAL trouvé. Plot non créé.")