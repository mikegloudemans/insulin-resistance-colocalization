#head -n -2 /users/bballiu/MoMeIR/data/IR_GWAS_Sumstats/ISI_GENESIS_adjAgeSex_JCI_paper_MetaAnalysis_Europeans_hg19.txt > /users/mgloud/projects/insulin_resistance/data/gwas/raw/sumstats/genesis_with_alleles/ISI_GENESIS_adjAgeSex_JCI_paper_MetaAnalysis_Europeans_hg19.txt
#head -n -2 /users/bballiu/MoMeIR/data/IR_GWAS_Sumstats/ISI_GENESIS_adjAgeSexBMI_JCI_paper_MetaAnalysis_Europeans_hg19.txt > /users/mgloud/projects/insulin_resistance/data/gwas/raw/sumstats/genesis_with_alleles/ISI_GENESIS_adjAgeSexBMI_JCI_paper_MetaAnalysis_Europeans_hg19.txt

#python /users/mgloud/projects/gwas-download/munge/custom_munge.py /users/mgloud/projects/insulin_resistance/scripts/format_gwas/munge_all_hg19.munge.config
python /users/mgloud/projects/gwas-download/munge/custom_munge.py /users/mgloud/projects/insulin_resistance/scripts/pre_coloc/format_gwas/munge_all_hg38.munge.config

