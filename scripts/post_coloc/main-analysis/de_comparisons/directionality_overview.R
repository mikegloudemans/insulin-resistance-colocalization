require(dplyr)

#data = read.table("output/post_coloc/de_genes/perturbation_by_coloc.txt", header=TRUE)
#data = read.table("output/post_coloc/de_genes/perturbation_by_single_coloc.txt", header=TRUE)

lower_is_better = c("data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz",
		    "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz")

higher_is_better = c("data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz",
		     "data/gwas/formatted/sumstats/hg38/MAGIC_ISI_Model_2_AgeSexBMI.txt/MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz")

# This is kind of hackish because it's doing a lot of double-counting, but still
# probably worth a look regardless
no_pert_data = data[,which(!(colnames(data) %in% c("pert_tissue", "perturbation", "pert_direction", "perturbation_increases_risk")))]
no_pert_data = no_pert_data[!duplicated(no_pert_data),]

direction_counts =  no_pert_data %>% group_by(hgnc) %>% summarize(num_high_causes_ir = sum((gwas_file %in% lower_is_better & higher_expression_higher_risk=="True") | (gwas_file %in% higher_is_better & higher_expression_higher_risk=="False")), num_low_causes_ir = sum((gwas_file %in% lower_is_better & higher_expression_higher_risk=="False") | (gwas_file %in% higher_is_better & higher_expression_higher_risk=="True")))

# Sometimes the direction is inconsistent across GWAS, but not that often
# and not in any way that indicates issues with the underlying intermediate data
sum((direction_counts$num_high_causes_ir > 0) & (direction_counts$num_low_causes_ir > 0))
inconsistent_genes = direction_counts[(direction_counts$num_high_causes_ir > 0) & (direction_counts$num_low_causes_ir > 0),]

# If we assign directions to every gene...we see no tendency for greater expression
# to lean towards either deleterious or protective phenotypes.
sum((direction_counts$num_high_causes_ir == 0) & (direction_counts$num_low_causes_ir > 0))
sum((direction_counts$num_high_causes_ir > 0) & (direction_counts$num_low_causes_ir == 0))

