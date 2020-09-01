# For part of figure 1A

require(dplyr)
d = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_snp_gene_pairs_prefiltering.txt", header=TRUE)
d$source_trait = sapply(d$source_trait, function(x) {s=strsplit(as.character(x), "/")[[1]]; s[length(s)]})

top_traits = d %>% group_by(ensembl) %>% summarize(IR = sum(source_trait %in% c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz")) > 0, FastInsu = !IR && sum(source_trait %in% c("FastInsu_adjBMI_MAGIC_Europeans.txt.gz") > 0), FastGlu = (sum(IR + FastInsu) == 0) && sum(source_trait %in% c("FastGlu_MAGIC_Europeans.txt.gz") > 0), T2D = (sum(IR + FastInsu + FastGlu) == 0) && sum(source_trait %in% c("Type-2-Diabetes_Spracklen_2020.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "T2D_Mahajan_Europeans.txt.gz") > 0), WHR = (sum(IR + FastInsu + FastGlu + T2D) == 0) && sum(source_trait %in% c("WHR-adj-BMI_GIANT_2018.txt.gz") > 0), TG = (sum(IR + FastInsu + FastGlu + T2D + WHR) == 0) && sum(source_trait %in% c("TG_GLGC_Expanded.txt.gz") > 0), HDL = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG) == 0) && sum(source_trait %in% c("HDL_GLGC_Expanded.txt.gz") > 0), BMI = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG + HDL) == 0) && sum(source_trait %in% c("BMI_GIANT_2018.txt.gz") > 0))

num_per_trait = colSums(top_traits[,-1])
sum(num_per_trait) # should equal total num genes with overlap
# it does for now at least


######################################################

coloc_tests = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", header=TRUE, fill=TRUE, sep="\t")

gene_tissues = coloc_tests %>% group_by(ensembl) %>% summarize(SubQ = sum(!(eqtl_short %in% c("Subcutaneous-Adipose-eQTL", "Subcutaneous-Adipose-sQTL"))) == 0, Visc = sum(!(eqtl_short %in% c("Visceral-Adipose-eQTL", "Visceral-Adipose-sQTL"))) == 0, Panc = sum(!(eqtl_short %in% c("Pancreas-eQTL", "Pancreas-sQTL"))) == 0, Liver = sum(!(eqtl_short %in% c("Liver-eQTL", "Liver-sQTL"))) == 0, Musc = sum(!(eqtl_short %in% c("Skeletal-Muscle-eQTL", "Skeletal-Muscle-sQTL"))) == 0, Shared=sum(SubQ + Visc + Panc + Liver + Musc) == 0)


tissue_nums = colSums(gene_tissues[,-1])
tissue_percs = colSums(gene_tissues[,-1]) / sum(colSums(gene_tissues[,-1])) * 100

sum(tissue_nums) # Should equal number of genes with QTL / GWAS overlap

#############################################################

single_coloc_at_locus = coloc_tests[coloc_tests$step1 %in% c("loci_0.1", "loci_2"),]

############################################################
d = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_gwas_snps.txt", header=TRUE)
d$source_trait = sapply(d$source_trait, function(x) {s=strsplit(as.character(x), "/")[[1]]; s[length(s)]})

locus_top_traits = d %>% group_by(locus) %>% summarize(IR = sum(source_trait %in% c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz")) > 0, FastInsu = !IR && sum(source_trait %in% c("FastInsu_adjBMI_MAGIC_Europeans.txt.gz") > 0), FastGlu = (sum(IR + FastInsu) == 0) && sum(source_trait %in% c("FastGlu_MAGIC_Europeans.txt.gz") > 0), T2D = (sum(IR + FastInsu + FastGlu) == 0) && sum(source_trait %in% c("Type-2-Diabetes_Spracklen_2020.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "T2D_Mahajan_Europeans.txt.gz") > 0), WHR = (sum(IR + FastInsu + FastGlu + T2D) == 0) && sum(source_trait %in% c("WHR-adj-BMI_GIANT_2018.txt.gz") > 0), TG = (sum(IR + FastInsu + FastGlu + T2D + WHR) == 0) && sum(source_trait %in% c("TG_GLGC_Expanded.txt.gz") > 0), HDL = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG) == 0) && sum(source_trait %in% c("HDL_GLGC_Expanded.txt.gz") > 0), BMI = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG + HDL) == 0) && sum(source_trait %in% c("BMI_GIANT_2018.txt.gz") > 0))

num_per_trait = colSums(locus_top_traits[,-1])
trait_percs = colSums(locus_top_traits[,-1]) / sum(colSums(locus_top_traits[,-1])) * 100
sum(num_per_trait) # should equal total num loci with overlap


###########################################################################################

locus_tissues = coloc_tests %>% group_by(locus) %>% summarize(SubQ = sum(!(eqtl_short %in% c("Subcutaneous-Adipose-eQTL", "Subcutaneous-Adipose-sQTL"))) == 0, Visc = sum(!(eqtl_short %in% c("Visceral-Adipose-eQTL", "Visceral-Adipose-sQTL"))) == 0, Panc = sum(!(eqtl_short %in% c("Pancreas-eQTL", "Pancreas-sQTL"))) == 0, Liver = sum(!(eqtl_short %in% c("Liver-eQTL", "Liver-sQTL"))) == 0, Musc = sum(!(eqtl_short %in% c("Skeletal-Muscle-eQTL", "Skeletal-Muscle-sQTL"))) == 0, Shared=sum(SubQ + Visc + Panc + Liver + Musc) == 0)
