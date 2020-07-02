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

g_a = NULL

##########################################################
# Bar chart by number candidates, number coloc
##########################################################

# Probably best to map colors from part A onto colors from part B

number_subcats = coloc_res %>% group_by(step1) %>% summarize(count = length(unique(locus)))

number_subcats$step1 = factor(x=number_subcats$step1, levels = c("loci_0.1", "loci_0.2", "loci_1",
						   "loci_2", "loci_3", "loci_4.1",
						   "loci_4.2", "loci_5"),
			labels = c("one candidate, one strong coloc", "one candidate, one weak coloc",
				   "one candidate, no coloc", "multi candidate, one strong coloc",
				   "multi candidate, multi strong coloc", "multi candidate, one weak coloc",
				   "multi candidate, multi weak coloc", "multi candidate, no coloc"))

ordering = rev(c(5,4,7,6,8,1,2,3))
category_color_scheme_ordered = category_color_scheme[ordering]
number_subcats$step1 = factor(x=number_subcats$step1, levels=number_subcats$step1[ordering])

# ^ This should probably be done with a multi-leveled legened on the axis -- one for num candidates,
# one for num genes

g_b <- ggplot(data=number_subcats, aes(x=step1, y=count, fill=step1)) +
	geom_bar(stat="identity", color="black") +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	scale_fill_manual(values = category_color_scheme_ordered) +
	geom_vline(xintercept=3.5, linetype="dashed") +
	theme(legend.position = "none") +
	theme(axis.title.y = element_blank())


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
number_subcats$step2 = factor(number_subcats$step2, levels = rev(c(nonpancreas_order, "Other", "Pancreas")))


g_c <- ggplot(data=number_subcats, aes(x=step2, y=count, fill=pancreas)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	scale_fill_manual(values=pancreas_presence) +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	theme(legend.position = "none") +
	theme(axis.title.y = element_blank())




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
			labels = rev(c("Insulin Resistance", "FastGlu, FastIns, T2D",
				   "Waist-Hip-Ratio", "HDL, TG",
				   "BMI")))

g_d <- ggplot(data=number_subcats, aes(x=step3, y=count, fill=pancreas)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	scale_fill_manual(values=pancreas_presence) +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	theme(axis.title.y = element_blank())


##########################################################
# (Maybe stratified by VEP)
##########################################################

##########################################################
# Now put them all together
##########################################################

g = plot_grid(g_a, g_b, g_c, g_d, labels = c('A', 'B', 'C', 'D'))

ggsave("output/post_coloc/plots/figure2/figure2.pdf", height=8, width=8)
