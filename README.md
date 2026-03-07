# BILL Project ( Bioinformatics Learning Lab )
## Introduction
This repository is dedicated to the BILL project of Group 5, conducted at the University of Montpellier. 
This project brings together students from the Bioinformatics and IMHE master’s programs. The data used are VCF files derived from viral cultures of Cyprinid herpesvirus 3 (KHV-U strain).


## USAGE

### Usage: DV_distribution_plot_sv scripts / QUAL_distribution_plot_snp (Python)
DV or QUAL distribution scripts are run from the Python folder (cd Python) using the following command:

```bash
python3 <script_name.py> <path_to_raw_VCF_folder>
```
The plots can be viewed in the Plots folder and enabled setting thresholds for the subsequent scripts.


### Usage: Extract_variant_filtered scripts (Python)
Filtered variant extraction scripts are run from the Python folder (cd Python) using the following command:

```bash
python3 <script_name.py> <path_to_raw_VCF_folder>
```

Output : All the filtered files are now located in <Files/vcf_filtered_sv> or <Files/vcf_filtered_snp>

### Usage: build_summary_table_sv_with_ORF.py
This script builds a summary table of structural variants (SV) from filtered VCF files (located in the vcf_filtered_sv folder), annotates each variant by ORF using the CyHV-3 GFF3 genome annotation, and export CSV file.

From the Python directory, use the command:
```bash
python3 <build_summary_table_sv_with_ORF.py>
```

Output : A detailed CSV of all variants with ORF annotation will be created in:
Files/plot_files/summary_table_sv_with_ORF.csv
This CSV can be used for downstream plotting and analysis.
