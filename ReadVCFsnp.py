#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
import os


argv=sys.argv
if len(argv)<2: #si l'argument a une taille supérieure à 2
    print("Usage: python ReadVCFsv.py <input.vcf>")
    sys.exit(1) 

vcf_folder=argv[1]
List_QUAL=[]
for filename in os.listdir(vcf_folder):
         if filename.endswith(".snp.vcf"):
            filepath = os.path.join(vcf_folder, filename)
            with open(filepath) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields=line.strip().split("\t")
                    chrom=fields[0]
                    pos=fields[1]
                    id=fields[2]
                    qual=fields[5]
                    qual = float(fields[5])
                    List_QUAL.append(qual)
print(len(List_QUAL))
print(f"QUAL: {List_QUAL}")

               # if qual >= 10:
                #    print(f"Chrom: {chrom}, Pos: {pos}, QUAL: {qual} → Haute qualité")
                 #   var_id = f"{chrom}_{pos}"
                  #  print(var_id)
output_file = "QUAL_list.csv"
with open(output_file, "w") as out:
    out.write(",".join(str(q) for q in List_QUAL))

print(f"Liste QUAL enregistrée dans {output_file} (valeurs séparées par des virgules)")