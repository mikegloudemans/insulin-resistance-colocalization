# For part of figure 1A

require(dplyr)
d = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_snp_gene_pairs_prefiltering.txt", header=TRUE)

d$source_trait = sapply(d$source_trait, function(x) {s=strsplit(as.character(x), "/")[[1]]; s[length(s)]})

top_traits = d %>% group_by(ensembl) %>% summarize(IR = sum(source_trait %in% c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz")) > 0, FastInsu = !IR && sum(source_trait %in% c("FastInsu_adjBMI_MAGIC_Europeans.txt.gz") > 0), FastGlu = (sum(IR + FastInsu) == 0) && sum(source_trait %in% c("FastGlu_MAGIC_Europeans.txt.gz") > 0), T2D = (sum(IR + FastInsu + FastGlu) == 0) && sum(source_trait %in% c("Type-2-Diabetes_Spracklen_2020.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "T2D_Mahajan_Europeans.txt.gz") > 0), WHR = (sum(IR + FastInsu + FastGlu + T2D) == 0) && sum(source_trait %in% c("WHR-adj-BMI_GIANT_2018.txt.gz") > 0), TG = (sum(IR + FastInsu + FastGlu + T2D + WHR) == 0) && sum(source_trait %in% c("TG_GLGC_Expanded.txt.gz") > 0), HDL = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG) == 0) && sum(source_trait %in% c("HDL_GLGC_Expanded.txt.gz") > 0), BMI = (sum(IR + FastInsu + FastGlu + T2D + WHR + TG + HDL) == 0) && sum(source_trait %in% c("BMI_GIANT_2018.txt.gz") > 0))

num_per_trait = colSums(top_traits[,-1])
sum(num_per_trait) # should equal total num genes with overlap
# it does for now at least

