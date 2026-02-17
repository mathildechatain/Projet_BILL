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
List_SUPPORT=[]
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
        SUPPORT=info2[4].split("=")[1]
        List_VAF.append(VAF)
        List_SUPPORT.append(SUPPORT)

    print(f"VAF: {List_VAF}")
    print(f"SUPPORT : {List_SUPPORT}")
