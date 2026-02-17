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
        qual = float(fields[5])
        List_QUAL.append(qual)

        if qual >= 10:
            print(f"Chrom: {chrom}, Pos: {pos}, QUAL: {qual} → Haute qualité")
            var_id = f"{chrom}_{pos}"
            print(var_id)