#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np


argv=sys.argv
if len(argv)<2: #si l'argument a une taille supérieure à 2
    print("Usage: python ReadVCFsv.py <input.vcf>")
    sys.exit(1) 

vcf_file=argv[1]
List_QUAL=[]
List_pos=[]
with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields=line.strip().split("\t")
        chrom=fields[0]
        pos=fields[1]
        id=fields[2]
        qual=fields[5]
        print(f"Chrom: {chrom}, Pos: {pos}, ID: {id}, QUAL: {qual}") 
        List_QUAL.append(float(qual)) 
        List_pos.append(pos)
    print(f"QUAL: {List_QUAL}")

fig, ax = plt.subplots()  # Il manquait les parenthèses
ax.plot(List_pos,List_QUAL)  # ajout de marqueurs pour mieux visualiser
ax.set_xlabel("position des variants")
ax.set_ylabel("QUAL")

ax.set_title("QUAL des variants dans le VCF")
plt.show()  #