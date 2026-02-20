import os
import glob
from pathlib import Path
import sys

# Check if the user provided the folder as an argument
if len(sys.argv) < 2:
    print("Usage: python filter_snp_vcf.py <dossier_vcf>") # else print an error message and exit
    sys.exit(1)

#the folder containing the VCF files is passed as an argument
vcf_folder = sys.argv[1]  
# Create the output folder, named vcf_filtered, if it doesn't exist
output_folder = os.path.join(vcf_folder, "vcf_filtered") 
Path(output_folder).mkdir(parents=True,exist_ok=True) 

# Function to filter the VCF file based on the DV (index 9, 4th value of FORMAT)
def filter_SV_vcf(input_file, output_file):

    # Open the input VCF file for reading and the output file for writing
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # read the input file line by line
        for line in infile:
            if line.startswith("#"):
                # keep header lines unchanged in the output file (columns names)
                outfile.write(line)
            else:
                # Split the line into columns using tab as a delimiter
                cols = line.strip().split("\t")
                try:
                    # Acces to the DV column (index 9, 4th value of FORMAT) 
                    DV = float(cols[9].split(":")[3]) 
                    # Filter the line based on the DV value (DV>=200)
                    if DV >= 200:
                        # IF the condition is met, write the line to the output file
                        outfile.write(line)
                except ValueError:
                    continue

# Find all VCF files in the input folder that end with "sniffles.vcf", that correspond to the SV variants
vcf_files = glob.glob(os.path.join(vcf_folder, "*sniffles.vcf"))

# Iterate over each VCF file and apply the filtering function
for vcf_file in vcf_files:
    # Generate the output file name by adding "_SV_filtered" to the base name of the input file
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    # Create the full path for the output file
    output_file = os.path.join(output_folder, f"{base_name}_SV_filtered.vcf")
    # Call the filtering function for each VCF file
    filter_SV_vcf(vcf_file, output_file)

# Print a message to indicate that the filtering process is completed
print("SV filtering completed!")

