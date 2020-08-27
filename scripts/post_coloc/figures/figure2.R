require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

source("scripts/post_coloc/figures/color_scheme.R")

# Plotting key info 

coloc_file = 'output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt'

coloc_res = read.table(coloc_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
coloc_res$base_gwas_file = sapply(coloc_res$base_gwas_file, function(x) {s = unlist(as.vector(strsplit(x, "/"))); s[length(s)]})

###########################################################
# Figure showing different colocalization types for "step
# 1": number candidate genes, number coloc genes
###########################################################

# old data -- I now will do this manually

# Probably best to map colors from part A onto colors from part B

number_subcats = coloc_res %>% group_by(step1) %>% summarize(count = length(unique(locus)))

number_subcats$step1 = factor(x=number_subcats$step1, levels = c("loci_0.1", "loci_0.2", "loci_1",
						   "loci_2", "loci_3", "loci_4.1",
						   "loci_4.2", "loci_5"),
			labels = c("one candidate, one coloc", "one candidate, one weak coloc",
				   "one candidate, no coloc", "multi candidate, one coloc",
				   "multi candidate, multi coloc", "multi candidate, one weak coloc",
				   "multi candidate, multi weak coloc", "multi candidate, no coloc"))

ordering = rev(c(5,4,7,6,8,1,2,3))
category_color_scheme_ordered = category_color_scheme[ordering]
number_subcats$step1 = factor(x=number_subcats$step1, levels=number_subcats$step1[ordering])

# ^ This should probably be done with a multi-leveled legened on the axis -- one for num candidates,
# one for num genes

g_a <- ggplot(data=number_subcats, aes(x=step1, y=count, fill=step1)) +
	geom_bar(stat="identity", color="black") +
	geom_text(color="black", aes(label=count, y=count+16), size=3) +
	coord_flip() +
	ylab("# GWAS loci") +
	theme_minimal() +
	scale_fill_manual(values = category_color_scheme_ordered) +
	geom_vline(xintercept=3.5, linetype="dashed") +
	theme(legend.position = "none") +
	theme(axis.title.y = element_blank())

##########################################################
# Bar chart by number candidates, number coloc
##########################################################

##########################################################
# Bar chart by number colocs in each tissue type
##########################################################

coloc_res_copy = coloc_res

# For now, simplify it so it's either pancreas or not pancreas
coloc_res_copy$step2[(coloc_res_copy$step2 %in% c("AdpS", "AdpV", "AdpV+S", "Liver", "Musk")) & (coloc_res_copy$pancreas == "weak")] = "Other"
coloc_res_copy$pancreas[coloc_res_copy$pancreas != "none"] = "present"
coloc_res_copy$pancreas[coloc_res_copy$pancreas == "none"] = "absent"
coloc_res_copy$pancreas = factor(coloc_res_copy$pancreas, levels=c("present", "absent"))

number_subcats = coloc_res_copy %>% group_by(step2, pancreas) %>% summarize(count = length(unique(locus)))
number_subcats = number_subcats[number_subcats$step2 != "None",]

nonpancreas_order = number_subcats$step2[order(number_subcats$count[number_subcats$step2 %in% c("AdpS", "AdpV", "AdpV+S", "Liver", "Musk")])]
tissue_types = c("Adipose (Subcutaneous)", "Adipose (Visceral)", "Adipose (Both)", "Liver", "Muscle (Skeletal)", "Shared, Pancreas Absent", "Shared, Pancreas Present", "Pancreas")
number_subcats$step2 = factor(tissue_types, levels = rev(tissue_types))

tissue_colors = rev(c("#FF6600", # Adipose-Subcutaneous
		  "#FFAA00", # Adipose-Visceral
		  "#d45500", # Adipose-Both
		  "#AABB66", # Liver
		  "#AAAAFF", # Skeletal-Muscle
		  "#c2c2c2", # Shared, Pancreas Absent
		  "#AF9079", # Shared, Pancreas Present
		  "#995522" # Pancreas
		 ))

g_b <- ggplot(data=number_subcats, aes(x=step2, y=count, fill=step2)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	geom_text(aes(label = count, y = count + 10))+
	scale_fill_manual(values=tissue_colors) +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	theme(legend.position = "none") +
	theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())




##########################################################
# Bar chart by number colocs in different coloc tiers
##########################################################

coloc_res_copy = coloc_res

# For now, simplify it so it's either pancreas or not pancreas
coloc_res_copy$pancreas[coloc_res_copy$pancreas != "none"] = "present"
coloc_res_copy$pancreas[coloc_res_copy$pancreas == "none"] = "absent"
coloc_res_copy$pancreas = factor(coloc_res_copy$pancreas, levels=c("present", "absent"))

number_subcats = coloc_res_copy %>% group_by(step3, pancreas) %>% summarize(count = length(unique(locus)))
number_subcats = number_subcats[number_subcats$step3 != "None",]


number_subcats$step3 = factor(x=number_subcats$step3, levels = rev(c("Tier_1", "Tier_2", "Tier_3",
						   "Tier_4", "Tier_5")),
			labels = rev(c("Insulin resistance", "Fasting glucose\nFasting insulin\nT2D",
				   "Waist-hip ratio", "HDL\nTriglycerides",
				   "BMI")))

g_c <- ggplot(data=number_subcats, aes(x=step3, y=count, fill=pancreas)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	geom_text(aes(label = count, y = count + 6, x = as.numeric(step3) + ifelse(pancreas=="present", -0.2, 0.2)))+
	scale_fill_manual(values=c("grey", "white"),
			  name = "pancreas\ncolocalization",
			  guide = guide_legend(reverse = TRUE)) +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


##########################################################
# (Maybe stratified by VEP)
##########################################################

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

##########################################################
# Now put them all together
##########################################################

g = plot_grid(g_a, g_b, g_c, g_d, labels = c('A', 'B', 'C', 'D'))

ggsave("output/post_coloc/plots/figure2/figure2.pdf", height=10, width=10)
