import os
import glob
from pathlib import Path
import sys

# check if the user provided the folder as an argument
if len(sys.argv) < 2:
    print("Usage: python filter_snp_vcf.py <dossier_vcf>") #else print an error message and exit    
    sys.exit(1)

#the folder containing the VCF files is passed as an argument
vcf_folder = sys.argv[1]  
# Create the output folder, named vcf_filtered, if it doesn't exist
output_folder = os.path.join(vcf_folder, "vcf_filtered")
Path(output_folder).mkdir(exist_ok=True)

# Function to filter SNPs based on our criterion QUAL > 10, which takes as input a VCF file and an output file
def filter_snp_vcf(input_file, output_file):

    # Open the input VCF file and the output file to write the filtered results
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Iterate over each line in the input file
        for line in infile:
            if line.startswith("#"):
                # Keep all header lines (columns names) unchanged in the output file
                outfile.write(line)
            else:
                # Split the line into columns to access the QUAL column
                cols = line.strip().split("\t")
                try:
                    # Access the QUAL column (index 5) and check if it's more than 10
                    qual = float(cols[5])
                    if qual > 10:
                        # If QUAL is greater than 10, write the line to the output file
                        outfile.write(line)
                except ValueError:
                    continue

# Processing all VCF files
# Use glob to find all VCF files in the specified folder that end with "snp.vcf", which correspond to the SNP variants
vcf_files = glob.glob(os.path.join(vcf_folder, "*snp.vcf"))
# Iterate over each VCF file found and apply the filtering function
for vcf_file in vcf_files:
    # Generate the output file name by adding "_SNP_filtered" to the base name of the input file
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    # Create the full path for the output file
    output_file = os.path.join(output_folder, f"{base_name}_SNP_filtered.vcf")
    # Call the filtering function for each VCF file
    filter_snp_vcf(vcf_file, output_file)

# Print a message to indicate that the filtering process is completed
print("SNP filtering completed!")

