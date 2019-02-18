# Sort dbSNP so we can easily join it with other files
zcat /mnt/lab_data/montgomery/shared/dbSNP/hg19_snp150.txt.gz | cut -f2,4,5 | sort -k3,3 | gzip > /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz

# Format individual GWAS results for use with coloc pipeline
bash prep_2hrGlu_MAGIC_Europeans.sh
bash prep_BMI_GIANT.sh
bash prep_CHD_CARDIoGRAMplusC4D_Mixed.sh
bash prep_FastGlu_MAGIC_Europeans.sh
bash prep_FastInsu_adjBMI_MAGIC_Europeans.sh
bash prep_FastInsu_MAGIC_Europeans.sh
bash prep_IS_METASTROKE_Mixed.sh
bash prep_MI_adjBMI_GENESISGuardian_Mixed.sh
bash prep_MI_GENESISGuardian_Mixed.sh
bash prep_T2D_DIAGRAM_European.sh
bash prep_WHRadjBMI_GIANT.sh
bash prep_ISI1_MAGIC.sh
bash prep_ISI2_MAGIC.sh


# Put into format for Bosh's web app
python format_for_bosh.py
