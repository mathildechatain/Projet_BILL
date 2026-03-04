#!/bin/bash

# Créer le dossier de sortie si nécessaire
mkdir -p ../Files/plot_files

# Fichier de sortie
output="../Files/plot_files/count_variants_sv.tsv"

# Écrire l'en-tête
echo -e "File\tN_variants" > "$output"

# Boucle sur tous les fichiers VCF merge
for f in ../Files/vcf_merged_sv/merged_*.vcf; do
  # Nom simplifié sans "merged_" et remplacer le point par "_"
  base=$(basename "$f" .vcf)
  file_name=$(echo "$base" | sed 's/^merged_//' | sed 's/\./_/')  # ex: P25_chaud

  # Trouver la ligne "#CHROM"
  chrom_line=$(grep -n "^#CHROM" "$f" | cut -d: -f1)

  # Compter le nombre de variants après #CHROM
  n_variants=$(tail -n +$((chrom_line+1)) "$f" | wc -l)

  # Ajouter au fichier de sortie
  echo -e "${file_name}\t${n_variants}" >> "$output"
done

echo "Résumé écrit dans $output"