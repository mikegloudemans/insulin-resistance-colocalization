require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

##########################################################
# (Maybe stratified by VEP)
##########################################################

# This is now part of the supplement, if that

annotations = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_1kg_filtered_2020-05-11.txt"
coloc_categories = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_classification_2020-05-11.txt"

summary_out_table = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/locus_vep_consequences.txt"

# for comparison
# pre_filtered_annotations = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_2020-05-11.txt"

effects = read.table(annotations)
colnames(effects)[c(3,4)] = c("locus", "effect")


effects$effect = factor(effects$effect,
				    levels = c("3PRIME_UTR", "5PRIME_UTR", "CANONICAL_SPLICE", "DOWNSTREAM", "INTERGENIC", "INTRONIC", 
					       "NONCODING_CHANGE", "NON_SYNONYMOUS", "REGULATORY", "SPLICE_SITE", "STOP_GAINED", 
					       "SYNONYMOUS", "UPSTREAM"),
				    labels = c("3' UTR", "5' UTR", "inside splice motif", "downstream", "intergenic", "intronic", "noncoding transcript", "non-synonymous", "regulatory element", "near splice motif", "stop-gain", "synonymous", "upstream"))

effects$impact = ""
effects$impact[effects$effect %in% c("inside splice motif", "stop-gain")] = "high"
effects$impact[effects$effect %in% c("non-synonymous")] = "moderate"
effects$impact[effects$effect %in% c("near splice motif", "synonymous")] = "low"
effects$impact[effects$effect %in% c("3' UTR", "5' UTR", "downstream", "intergenic", "intronic", "noncoding transcript", "regulatory element","upstream")] = "modifier"
effects$impact = factor(effects$impact, levels = c("high", "moderate", "low", "modifier"))

coloc_data = read.table(coloc_categories, header=TRUE, sep="\t")
coloc_only = coloc_data %>% filter(step2 != "None")
coloc_effects = effects %>% filter(locus %in% unique(coloc_only$locus))

coloc_effect_counts = coloc_effects %>% group_by(effect) %>% summarize(count=length(unique(locus)), percent = count/length(unique(coloc_effects$locus))*100, impact=impact[1])
coloc_effect_counts = coloc_effect_counts %>% arrange(percent)
coloc_effect_counts$effect = factor(coloc_effect_counts$effect, levels = coloc_effect_counts$effect)

g_d1 <- ggplot(data=coloc_effect_counts, aes(x=effect, y=percent, fill=impact)) +
	geom_bar(stat="identity", color="black", position="dodge") +
	geom_text(color="black", aes(label=count, y=percent+6), size=3) +
	coord_flip() +
	theme_minimal() +
	scale_fill_manual(values = c("red", "yellow", "green", "grey85")) +
	#theme(legend.position = "none") +
	theme(axis.title.y = element_blank())

coloc_effect_impact_counts = coloc_effects %>% group_by(impact) %>% summarize(count=length(unique(locus)), percent = count/length(unique(coloc_effects$locus))*100)
coloc_effect_impact_counts = coloc_effect_impact_counts %>% filter(impact != "modifier")
coloc_effect_impact_counts = coloc_effect_impact_counts %>% arrange(percent)

# Also include one with aggregate info
g_d2_pre <- ggplot(data=coloc_effect_impact_counts, aes(x=impact, y=percent, fill=impact)) +
	geom_bar(stat="identity", color="black", position="dodge") +
	geom_text(color="black", aes(label=count, y=percent+6), size=3) +
	coord_flip() +
	ylab("% of coloc loci with annotation") +
	ylim(c(0,100)) +
	theme_minimal() +
	scale_fill_manual(values = c("red", "yellow", "green", "grey85")) +
	theme(legend.position = "none") +
	theme(axis.title.y = element_blank())

g_d2 = plot_grid(NULL, g_d2_pre, NULL, ncol=3, rel_widths=c(0.13,0.69,0.18))

g_d = plot_grid(g_d1, g_d2, nrow=2, rel_heights=c(2.2,1))

ggsave("output/post_coloc/plots/vep-consequences.pdf", height=4, width=15)
