## GIANT

# Lift BMI hg18 to hg19
echo -e "rsid\tchr\tsnp_pos\tgene\talt\tref\talt_freq\tbeta\tse\tvariance_explained\tN\tpvalue" \
	> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/BMI_GIANT_Europeans_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f4,5 --complement >> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tgene\talt\tref\talt_freq\tbeta\tse\tvariance_explained\tN\tpvalue" \
	> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/BMI_GIANT_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f4,5,15,16 --complement >> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Mixed_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\ttrait\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Europeans_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/WHRadjBMI_GIANT_Europeans_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f5,6 --complement >> /users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\ttrait\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/WHRadjBMI_GIANT_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f5,6 --complement >> /users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Mixed_Hits_hg19.txt


## MAGIC

echo -e "rsid\tchr\tsnp_pos\tcategory\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN\tBMIAdj_beta\tBMIAdj_SE\tBMIAdj_pvalue\tBMIAdj_N" \
	> /users/mgloud/projects/insulin_resistance/data/2hrGlu_MAGIC_Europeans_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/2hrGlu_MAGIC_Europeans_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5,6 --complement >> /users/mgloud/projects/insulin_resistance/data/2hrGlu_MAGIC_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tcategory\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN\tBMIAdj_beta\tBMIAdj_SE\tBMIAdj_pvalue\tBMIAdj_N" \
	> /users/mgloud/projects/insulin_resistance/data/FastGlu_MAGIC_Europeans_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/FastGlu_MAGIC_Europeans_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5,6 --complement >> /users/mgloud/projects/insulin_resistance/data/FastGlu_MAGIC_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tcategory\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN\tBMIAdj_beta\tBMIAdj_SE\tBMIAdj_pvalue\tBMIAdj_N" \
	> /users/mgloud/projects/insulin_resistance/data/FastInsu_MAGIC_Europeans_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/FastInsu_MAGIC_Europeans_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5-6 --complement >> /users/mgloud/projects/insulin_resistance/data/FastInsu_MAGIC_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tcategory\tgene\tnearer_gene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN_studies\tN_samples" \
	> /users/mgloud/projects/insulin_resistance/data/FastInsu_adjBMI_MAGIC_Europeans_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/FastInsu_adjBMI_MAGIC_Europeans_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5-6 --complement >> /users/mgloud/projects/insulin_resistance/data/FastInsu_adjBMI_MAGIC_Europeans_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tcategory\teffect_allele\tnon_effect_allele\tgene\tstatus\tsignals\tclassification\tpvalue" \
	> /users/mgloud/projects/insulin_resistance/data/HbA1c_MAGIC_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/HbA1c_MAGIC_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f4-5 --complement >> /users/mgloud/projects/insulin_resistance/data/HbA1c_MAGIC_Mixed_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tgene\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN\tnotes" \
	> /users/mgloud/projects/insulin_resistance/data/ISI2_MAGIC_Europeans_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/ISI2_MAGIC_Europeans_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/insulin_resistance/data/ISI2_MAGIC_Europeans_Hits_hg19.txt

# GCLC

echo -e "rsid\tchr\tsnp_pos\tgene\tlocus\ttrait\tmaf\talleles\ta1_effect\tjoint_N_in_thousands\tjoint_p_value" \
	> /users/mgloud/projects/insulin_resistance/data/HDL_GCLC_Mixed_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/HDL_GCLC_Mixed_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5 --complement >> /users/mgloud/projects/insulin_resistance/data/HDL_GCLC_Mixed_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\tgene\tlocus\ttrait\tmaf\talleles\ta1_effect\tjoint_N_in_thousands\tjoint_p_value" \
	> /users/mgloud/projects/insulin_resistance/data/TG_GCLC_Mixed_Hits_hg19.txt
join -1 3 -2 2 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/TG_GCLC_Mixed_Hits.txt | tail -n +2 | sort -k2,2) | sort -k2,2 -k3,3n | cut -f5 --complement >> /users/mgloud/projects/insulin_resistance/data/TG_GCLC_Mixed_Hits_hg19.txt

# GENESIS
echo -e "rsid\tchr\tsnp_pos\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/insulin_resistance/data/MI_adjBMI_GENESISGuardian_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(cat /users/mgloud/projects/insulin_resistance/data/tophits/MI_adjBMI_GENESISGuardian_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/insulin_resistance/data/MI_adjBMI_GENESISGuardian_Mixed_Hits_hg19.txt

echo -e "rsid\tchr\tsnp_pos\teffect_allele\tnon_effect_allele\tea_freq\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/insulin_resistance/data/MI_GENESISGuardian_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(cat /users/mgloud/projects/insulin_resistance/data/tophits/MI_GENESISGuardian_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/insulin_resistance/data/MI_GENESISGuardian_Mixed_Hits_hg19.txt

# DIAGRAM
echo -e "rsid\tchr\tsnp_pos\tgene\ta1\ta2\tfreq1\tor1\tci1\tpvalue1\tN1\tor2\tci2\tpvalue2\tN2\ttstat\thet_pvalue" \
	> /users/mgloud/projects/insulin_resistance/data/T2D_DIAGRAM_European_Hits_hg19.txt
join -1 3 -2 5 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/T2D_DIAGRAM_European_Hits.txt | tail -n +2 | sort -t$'\t' -k5,5) | sort -k2,2 -k3,3n | cut -f4-6 --complement >> /users/mgloud/projects/insulin_resistance/data/T2D_DIAGRAM_European_Hits_hg19.txt

# METASTROKE
echo -e "rsid\tchr\tsnp_pos\tgene\ttrait\teffect_allele\tea_freq_discovery\tpvalue\tor" \
	> /users/mgloud/projects/insulin_resistance/data/IS_METASTROKE_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/IS_METASTROKE_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/insulin_resistance/data/IS_METASTROKE_Mixed_Hits_hg19.txt

# CARDIOGRAM
echo -e "rsid\tchr\tsnp_pos\tgene\talleles\tea_freq\tor\tpvalue\tor_recessive\tpvalue_recessive\ttype" \
	> /users/mgloud/projects/insulin_resistance/data/CHD_CARDIoGRAMplusC4D_Mixed_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(sed s//\\n/g /users/mgloud/projects/insulin_resistance/data/tophits/CHD_CARDIoGRAMplusC4D_Mixed_Hits.txt | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | cut -f5 --complement >> /users/mgloud/projects/insulin_resistance/data/CHD_CARDIoGRAMplusC4D_Mixed_Hits_hg19.txt


