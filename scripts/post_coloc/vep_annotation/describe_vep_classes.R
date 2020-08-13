% with each annotation
% with / without a "serious" annotation

require(dplyr)
require(ggplot2)

annotations = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_1kg_filtered_2020-05-11.txt"
coloc_categories = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_classification_2020-05-11.txt"

summary_out_table = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/locus_vep_consequences.txt"

# for comparison
# pre_filtered_annotations = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_2020-05-11.txt"

effects = read.table(annotations)
colnames(effects)[c(3,4)] = c("locus", "effect")
effects$importance = ifelse(effects$effect %in% c("STOP_GAINED", "CANONICAL_SPLICE", "NON_SYNONYMOUS"), "serious", "minor")
coloc_data = read.table(coloc_categories, header=TRUE, sep="\t")
coloc_only = coloc_data %>% filter(step2 != "None")
coloc_effects = effects %>% filter(locus %in% unique(coloc_only$locus))

##########################################################
# Display percent of loci with each annotation
##########################################################

effect_counts = effects %>% group_by(effect) %>% summarize(count=length(unique(locus)), percent = count/length(unique(effects$locus))*100)
effect_importance_counts = effects %>% group_by(importance) %>% summarize(count=length(unique(locus)), percent = count/length(unique(effects$locus))*100)
coloc_effect_counts = coloc_effects %>% group_by(effect) %>% summarize(count=length(unique(locus)), percent = count/length(unique(coloc_effects$locus))*100)
coloc_effect_importance_counts = coloc_effects %>% group_by(importance) %>% summarize(count=length(unique(locus)), percent = count/length(unique(coloc_effects$locus))*100)
coloc_effect_counts = coloc_effect_counts %>% arrange(percent)
coloc_effect_importance_counts = coloc_effect_importance_counts %>% arrange(percent)

coloc_effect_counts$type = "coloc_loci"
coloc_effect_importance_counts$type = "coloc_loci"
effect_counts$type = "all_loci"
effect_importance_counts$type = "all_loci"

combined_effect_counts = rbind(coloc_effect_counts, effect_counts)
combined_effect_counts$effect = factor(combined_effect_counts$effect, levels = coloc_effect_counts$effect)
combined_effect_importance_counts = rbind(coloc_effect_importance_counts, effect_importance_counts)

g <- ggplot(data=combined_effect_counts, aes(x=effect, y=percent, fill=type)) +
	geom_bar(stat="identity", color="black", position="dodge") +
	geom_text(color="black", aes(label=count, y=percent+6, x=as.numeric(effect)+ifelse(as.character(type)=="all_loci", -0.22, 0.22)), size=3) +
	coord_flip() +
	ylab("% of loci with annotation") +
	theme_minimal() +
	#scale_fill_manual(values = category_color_scheme_ordered) +
	#theme(legend.position = "none") +
	theme(axis.title.y = element_blank())

g <- ggplot(data=combined_effect_importance_counts, aes(x=importance, y=percent, fill=type)) +
	geom_bar(stat="identity", color="black", position="dodge") +
	geom_text(color="black", aes(label=count, y=percent+6, x=as.numeric(importance)+ifelse(as.character(type)=="all_loci", -0.22, 0.22)), size=3) +
	coord_flip() +
	ylab("% of loci with annotation") +
	theme_minimal() +
	#scale_fill_manual(values = category_color_scheme_ordered) +
	#theme(legend.position = "none") +
	theme(axis.title.y = element_blank())


###############################################################
# Create simple table showing what annotations there are at 
# each locus
###############################################################

severe_effects = effects %>% group_by(locus) %>% summarize(coloc=locus[1] %in% coloc_only$locus, severe_consequences = gsub("^,|,$", "", paste(unique(ifelse(effect %in% c("NON_SYNONYMOUS", "STOP_GAINED", "CANONICAL_SPLICE"), as.character(effect), "")), collapse=",")))

write.table(severe_effects, file=summary_out_table, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

