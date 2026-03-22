#!/usr/bin/env bash

# option de sécurité
set -euo pipefail

# Récupère le chemin absolu du dossier racine du projet
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Fonction utilitaire pour exécuter une étape du pipeline
run_step() {
  echo
  echo "==> $1"   # Affiche le nom de l'étape
  shift        # Retire le premier argument (le message)
  (cd "$PROJECT_ROOT" && "$@")  # Exécute la commande à la racine du projet
}


# Génération du tableau récapitulatif des SV avec annotation ORF
run_step "Generate summary_table_sv_with_ORF CSV" \
  python3 src/Python/build_summary_table_sv_with_ORF.py

# Génération des comptages de variants structuraux (SV)
run_step "Generate SV counts CSV" \
  python3 src/Python/build_sv_sc_counts.py

# Visualisation des SV par échantillon et génération
run_step "Plot SV counts by sample and generation" \
  Rscript src/R/count_sv_by_echantillon_gen.R

# Génération du CSV pour la sélection des passages
run_step "Generate passage selection CSV" \
  python3 src/Python/build_csv_selection_passages.py

# Visualisation de la sélection des passages
run_step "Plot passage selection" \
  Rscript src/R/SV_selection_par_passage.R

# Analyse des variants partagés entre P25 et P27 + évolution de la VAF
run_step "Generate shared variants P25/P27 summary and plot" \
  python3 src/Python/shared_variants_P25_P27_vaf_evolution.py

# Génération des SV par ORF
run_step "Generate SV by ORF CSVs" \
  python3 src/Python/sv_by_orf_generation.py

# Visualisation des SV par ORF et génération
run_step "Plot SV by ORF and generation" \
  Rscript src/R/sv_by_orf_generation.R

# Visualisation globale des SV par génération
run_step "Plot SV per generation" \
  Rscript src/R/sv_by_generation.R

# Analyse du turnover des SV entre P25 et P27 (groupes chaud/froid)
run_step "Plot SV turnover between P25 and P27 for chaud/froid groups" \
  Rscript src/R/SV_turnover_P25_P27_chaud_froid.R

# Analyse des ORFs les plus impactés après choc thermique
run_step "Plot most impacted ORFs after thermal shock" \
  Rscript src/R/ORF_impact_choc_thermique.R

# Message final
echo
echo "pipeline SV terminé avec succès"