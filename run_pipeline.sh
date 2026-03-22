#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_step() {
  echo
  echo "==> $1"
  shift
  (cd "$PROJECT_ROOT" && "$@")
}

run_step "Generate summary_table_sv_with_ORF CSV" \
  python3 src/Python/build_summary_table_sv_with_ORF.py

run_step "Generate SV counts CSV" \
  python3 src/Python/build_sv_sc_counts.py

run_step "Plot SV counts by sample and generation" \
  Rscript src/R/count_sv_by_echantillon_gen.R

run_step "Generate passage selection CSV" \
  python3 src/Python/build_csv_selection_passages.py

run_step "Plot passage selection" \
  Rscript src/R/SV_selection_par_passage.R

run_step "Generate shared variants P25/P27 summary and plot" \
  python3 src/Python/shared_variants_P25_P27_vaf_evolution.py

run_step "Generate SV by ORF CSVs" \
  python3 src/Python/sv_by_orf_generation.py

run_step "Plot SV by ORF and generation" \
  Rscript src/R/sv_by_orf_generation.R

run_step "Plot SV per generation" \
  Rscript src/R/sv_by_generation.R

run_step "Plot SV turnover between P25 and P27 for chaud/froid groups" \
  Rscript src/R/SV_turnover_P25_P27_chaud_froid.R

run_step "Plot most impacted ORFs after thermal shock" \
  Rscript src/R/ORF_impact_choc_thermique.R

echo
echo "Pipeline completed successfully."
