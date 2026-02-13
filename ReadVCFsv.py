#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np


argv=sys.argv
if len(argv)<2: #si l'argument a une taille supérieure à 2
    print("Usage: python ReadVCFsv.py <input.vcf>")
    sys.exit(1) 

vcf_file=argv[1]
List_VAF=[]
List_pos=[]
with open(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields=line.strip().split("\t")
        chrom=fields[0]
        pos=fields[1]
        id=fields[2]
        filter=fields[6]
        info=fields[7]
        print(f"Chrom: {chrom}, Pos: {pos}, ID: {id}, Filter: {filter}, Info: {info}")    
        info2=info.split(";")
        VAF=info2[10].split("=")[1]
        List_VAF.append(VAF)
        List_pos.append(pos)
    print(f"VAF: {List_VAF}")

fig, ax = plt.subplots()  # Il manquait les parenthèses
ax.plot(List_pos,List_VAF)  # ajout de marqueurs pour mieux visualiser
ax.set_xlabel("position des variants")
ax.set_ylabel("VAF")
ax.set_title("VAF des variants dans le VCF")
plt.show()  #