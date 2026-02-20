#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import os

argv=sys.argv
if len(argv)<2: #si l'argument a une taille supérieure à 2
    print("Usage: python ReadVCFsv.py <input.vcf>")
    sys.exit(1) 

vcf_folder=argv[1]
List_DV=[]

for filename in os.listdir(vcf_folder):
         if filename.endswith("sniffles.vcf"):
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

output_file_support = "List_DV.txt"
with open(output_file_support, "w") as out:
    out.write(",".join(str(q) for q in List_DV))

print(f"Liste DV enregistrée dans {output_file_support} (valeurs séparées par des virgules)")
with open(output_file_support, "r") as f:
    line = f.readline()  # tout est sur une seule ligne
    DV_values = [float(x) for x in line.strip().split(",")]  # convertir en float

plt.figure(figsize=(10,6))
plt.hist(DV_values, bins=50, color='skyblue', edgecolor='black')
plt.title("Distribution des valeurs DV des variants")
plt.xlabel("DV ( Depth Variant)")
plt.ylabel("Nombre de variants")
plt.grid(axis='y', alpha=0.75)
plt.show()
