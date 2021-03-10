require(dplyr)

data = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", header=TRUE, sep="\t")

tissues = data %>% group_by(locus) %>% summarize(AdpS = sum(!(eqtl_file %in% c("data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz", "data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz")))==0,
				       AdpV = sum(!(eqtl_file %in% c("data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz", "data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz")))==0,
				       Liver = sum(!(eqtl_file %in% c("data/sqtls/gtex_v8/Liver.sQTLs.txt.gz", "data/eqtls/gtex_v8/Liver.allpairs.txt.gz.eQTLs.txt.gz")))==0,
				       Muscle = sum(!(eqtl_file %in% c("data/sqtls/gtex_v8/Muscle_Skeletal.sQTLs.txt.gz", "data/eqtls/gtex_v8/Muscle_Skeletal.allpairs.txt.gz.eQTLs.txt.gz")))==0,
				       Pancreas = sum(!(eqtl_file %in% c("data/sqtls/gtex_v8/Pancreas.sQTLs.txt.gz", "data/eqtls/gtex_v8/Pancreas.allpairs.txt.gz.eQTLs.txt.gz")))==0,
				       Shared = !(AdpS || AdpV || Liver || Muscle || Pancreas))
