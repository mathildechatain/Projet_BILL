# BILL Project ( Bioinformatics Learning Lab )
## Introduction
This repository is dedicated to the BILL project of Group 5, conducted at the University of Montpellier. 
This project brings together students from the Bioinformatics and IMHE master’s programs. The data used are VCF files derived from viral cultures of Cyprinid herpesvirus 3 (KHV-U strain).


## USAGE

### 1. Quality distribution plots
`DV_distribution_plot_sv` and `QUAL_distribution_plot_snp` scripts are run from `src/Python` using:

```bash
python3 <script_name.py> <path_to_raw_VCF_folder>
```

The plots are saved in `Plots/` and can be used to define filtering thresholds for downstream analyses.


### 2. Filtered variant extraction
Filtered variant extraction scripts are run from `src/Python` using:

```bash
python3 Extract_variant_filtered_SV.py <path_to_raw_VCF_folder>
python3 Extract_variant_filtered_SNP.py <path_to_raw_VCF_folder>
```

Output: filtered files are written in `Files/vcf_filtered_sv/` or `Files/vcf_filtered_snp/`.


### 3. Automated pipelines
After generating filtered VCFs, the full downstream analyses can be launched from the project root with:

```bash
./run_pipeline_SV.sh
./run_pipeline_snp.sh
```

`run_pipeline_SV.sh` generates all SV analysis CSV files and plots by executing the required Python and R scripts in the correct order.

`run_pipeline_snp.sh` generates all SNP analysis CSV files and plots by executing the required Python and R scripts in the correct order.
