
bosh_dir="/users/mgloud/projects/gwas/data/boshs_format"

#mv $bosh_dir/adhd.prepared.txt.bosh $bosh_dir/final/GWAS_ADHD_unknown_2000.txt
#mv $bosh_dir/adiponectin_adipogen.prepared.txt.bosh $bosh_dir/final/GWAS_Adiponectin_unknown_2000.txt
#mv $bosh_dir/age_at_menopause.prepared.txt.bosh $bosh_dir/final/GWAS_AgeAtMenopause_unknown_2000.txt
#mv $bosh_dir/aggression.prepared.txt.bosh $bosh_dir/final/GWAS_Aggression_unknown_2000.txt

mv $bosh_dir/asthma.prepared.txt.bosh $bosh_dir/final/GWAS_Asthma_Moffatt_2010.txt
mv $bosh_dir/rheumatoid_arthritis_Asian.prepared.txt.bosh $bosh_dir/final/GWAS_Rheumatoid-Arthritis-Asian_Okada_2014.txt
mv $bosh_dir/rheumatoid_arthritis_European.prepared.txt.bosh $bosh_dir/final/GWAS_Rheumatoid-Arthritis-European_Okada_2014.txt
mv $bosh_dir/rheumatoid_arthritis_Mixed.prepared.txt.bosh $bosh_dir/final/GWAS_Rheumatoid-Arthritis-Mixed_Okada_2014.txt


##
## Files processed as part of Bruna's coloc analysis
##
mv $bosh_dir/BMI_GIANT_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_BMI-Europeans_Locke_2015.txt
mv $bosh_dir/BMI_GIANT_Mixed_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_BMI-Mixed_Locke_2015.txt
mv $bosh_dir/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt.bosh $bosh_dir/final/GWAS_Coronary-Heart-Disease_Nikpay_2015.txt
mv $bosh_dir/IS_METASTROKE_Mixed.prepared.txt.bosh $bosh_dir/final/GWAS_Ischemic-Stroke_Malik_2016.txt
mv $bosh_dir/T2D_DIAGRAM_European.prepared.txt.bosh $bosh_dir/final/GWAS_Type-2-Diabetes_Scott_2017.txt
mv $bosh_dir/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_Waist-Hip-Ratio-adjusted-for-BMI-Europeans_Shungin_2015.txt
mv $bosh_dir/WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_Waist-Hip-Ratio-adjusted-for-BMI-Mixed_Shungin_2015.txt

# These ones don't have enough SNPs profiled
#mv $bosh_dir/2hrGlu_MAGIC_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_2hrGlu__2000.txt
#mv $bosh_dir/age_at_menarche.prepared.txt.bosh $bosh_dir/final/GWAS_AgeAtMenarche_unknown_2000.txt
#mv $bosh_dir/FastGlu_MAGIC_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_FastGlu_unknown_2000.txt
#mv $bosh_dir/FastInsu_MAGIC_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_FastInsu_unknown_2000.txt
#mv $bosh_dir/FastInsu_adjBMI_MAGIC_Europeans_AllSNPs.prepared.txt.bosh $bosh_dir/final/GWAS_FastInsu-adjBMI_unknown_2000.txt
#mv $bosh_dir/MI_adjBMI_GENESISGuardian_Mixed.prepared.txt.bosh $bosh_dir/final/GWAS_MI-adjBMI_unknown_2000.txt
#mv $bosh_dir/MI_GENESISGuardian_Mixed.prepared.txt.bosh $bosh_dir/final/GWAS_MI_unknown_2000.txt
