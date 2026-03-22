#!/usr/bin/env bash

# options de sécurité :
set -euo pipefail

# Récupère le chemin absolu du dossier racine du projet
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Fonction pour exécuter une étape Python
run_python_step() {
  echo
  echo "==> $1"   # Affiche le nom de l'étape
  shift        # Supprime le message
  (cd "$PROJECT_ROOT" && "$@")  # Exécute la commande à la racine du projet
}

# Fonction pour exécuter une étape R
run_r_step() {
  echo
  echo "==> $1"
  shift
  # Exécute les scripts R depuis le dossier src/R
  (cd "$PROJECT_ROOT/src/R" && Rscript "$@")
}



# Génération du tableau récapitulatif des SNP avec annotation ORF
run_python_step "Generate summary_table_snp_with_ORF CSV" \
  python3 src/Python/build_summary_table_snp_with_ORF.py

# Génération des comptages de SNP
run_python_step "Generate SNP counts CSV" \
  python3 src/Python/build_snp_sc_counts.py

# Visualisation des SNP par échantillon et génération
run_r_step "Plot SNP counts by sample and generation" \
  count_snp_by_repetition.R

# Génération du CSV pour la sélection des passages SNP
run_python_step "Generate SNP passage selection CSV" \
  python3 src/Python/build_snp_sc_selection_with_orf.py

# Visualisation de la sélection des passages SNP
run_r_step "Plot SNP passage selection" \
  SNP_selection_par_passage.R

# Visualisation des SNP par génération
run_r_step "Plot SNP per generation" \
  snp_by_generation.R

# Visualisation des SNP par ORF
run_r_step "Plot SNP by ORF" \
  SNP_par_ORF

# Message final
echo
echo "pipeline SNP terminé avec succès"