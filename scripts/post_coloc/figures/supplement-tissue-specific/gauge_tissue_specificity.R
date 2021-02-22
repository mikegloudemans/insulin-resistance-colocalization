require(dplyr)
require(reshape2)
require(ggplot2)

data = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", sep="\t", header=TRUE)

#tissue_specific = data %>% filter(!(step2 %in% c("None", "Other")))

tissue_diffs = function(cat_name, tissue_grep_query)
{
	adpv_specific = data %>% filter(step2 == cat_name)
	adpv_unique = adpv_specific %>% filter(clpp_mod >= 0.35)
	adpv_unique = unique(adpv_unique[c("locus", "gwas_short", "hgnc")])
	adpv_diffs = sapply(1:dim(adpv_unique)[1], function(x)
	       {
			dats = list()

			this_locus = adpv_unique[x,]

			sub_in = adpv_specific %>% filter(locus == this_locus$locus) %>% filter(gwas_short == this_locus$gwas_short) %>% filter(hgnc == this_locus$hgnc) %>% filter(grepl(tissue_grep_query, eqtl_short))
			sub_out = adpv_specific %>% filter(locus == this_locus$locus) %>% filter(gwas_short == this_locus$gwas_short) %>% filter(hgnc == this_locus$hgnc) %>% filter(!grepl(tissue_grep_query, eqtl_short))

			dats["target-tissue"] = mean(sub_in$clpp_mod)
			dats["other-tissue"] = mean(sub_out$clpp_mod)

			return(dats)
	       })
	adpv_diffs = adpv_diffs[,colSums(is.na(adpv_diffs)) == 0]
}

adpv_diffs = tissue_diffs("AdpV", "Visceral-Adipose")
adps_diffs =  tissue_diffs("AdpS", "Subcutaneous-Adipose")
adpsv_diffs =  tissue_diffs("AdpV+S", "Adipose")
pancreas_diffs =  tissue_diffs("Pancreas", "Pancreas")
musk_diffs =  tissue_diffs("Musk", "Muscle")
liver_diffs =  tissue_diffs("Liver", "Liver")

all_diffs = cbind(adpv_diffs, adps_diffs, adpsv_diffs, pancreas_diffs, musk_diffs, liver_diffs)
diff_mat = as.data.frame(list(target_tissue = unlist(all_diffs[1,]), other_tissues = unlist(all_diffs[2,])))

in_data = melt(diff_mat)
colnames(in_data) = c("tissue", "CLPP_mod")

g = ggplot(in_data, aes(x=CLPP_mod, fill=tissue)) +
       	geom_density(alpha=0.6) +
	scale_fill_discrete(name="", labels=c("coloc tissue", "other tested tissues"))+
	xlab("Colocalization probability % (CLPP-mod)")


ggsave("output/post_coloc/plots/supplement-tissue-specific/supplement-tissue-specific.pdf", width=6, height=4)
