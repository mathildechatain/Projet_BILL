mkdir -p ../Files/vcf_zipped_copy
# Créer le dossier de sortie si nécessaire

# Copier tous les fichiers .vcf depuis vcf_filtered vers le dossier de copie
cp ../Files/vcf_filtered_sv/*.vcf \
   ../Files/vcf_zipped_copy/
# Se placer dans le dossier de copie
cd ../Files/vcf_zipped_copy
# Compresser tous les fichiers .vcf
for file in *.vcf; do
    bgzip "$file"
    tabix -p vcf "${file}.gz"
done

mkdir -p ../Files/merged_sv_filtered

bcftools merge P25-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o ../Files/merged_sv_filtered/merged_P25.froid.vcf
bcftools merge P27-1.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-2.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-3.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-4.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-5.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o ../Files/merged_sv_filtered/merged_P27.froid.vcf
bcftools merge P25-6.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-7.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-8.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-9.trimed1000.sv_sniffles_SV_filtered.vcf.gz P25-10.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o ../Files/merged_sv_filtered/merged_P25.chaud.vcf
bcftools merge P27-6.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-7.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-8.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-9.trimed1000.sv_sniffles_SV_filtered.vcf.gz P27-10.trimed1000.sv_sniffles_SV_filtered.vcf.gz \
    -o ../Files/merged_sv_filtered/merged_P27.chaud.vcf