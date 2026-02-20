
mkdir -p ~/Bureau/Master/M1/S2/BILL/Ressources/vcf_zipped_copy

# Copier tous les fichiers .vcf depuis vcf_filtered vers le dossier de copie
cp ~/Bureau/Master/M1/S2/BILL/Ressources/vcf_filtered/*.vcf \
   ~/Bureau/Master/M1/S2/BILL/Ressources/vcf_zipped_copy/

# Se placer dans le dossier de copie
cd ~/Bureau/Master/M1/S2/BILL/Ressources/vcf_zipped_copy

# Compresser tous les fichiers .vcf
for file in *.vcf; do
    bgzip "$file"
    tabix -p vcf "${file}.gz"
done

bcftools merge P25-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o merged_P25.froid.vcf
bcftools merge P27-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o merged_P27.froid.vcf
bcftools merge P25-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    P27-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o merged_P25-27.froid.vcf
