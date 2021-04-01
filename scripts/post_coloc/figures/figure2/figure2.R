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

# This part was done manually, using the per-type counts given
# in the "classify results" section of the main analysis

##########################################################
# Bar chart by number colocs in each tissue type
##########################################################

coloc_res_copy = coloc_res

# For now, simplify it so it's either pancreas or not pancreas
coloc_res_copy$step2[(coloc_res_copy$step2 %in% c("AdpS", "AdpV", "AdpV+S", "Liver", "Musk")) & (coloc_res_copy$pancreas == "weak")] = "Other"
coloc_res_copy$pancreas[coloc_res_copy$pancreas != "none"] = "present"
coloc_res_copy$pancreas[coloc_res_copy$pancreas == "none"] = "absent"
coloc_res_copy$pancreas = factor(coloc_res_copy$pancreas, levels=c("absent", "present"))

number_subcats = coloc_res_copy %>% group_by(step2, pancreas) %>% summarize(count = length(unique(locus)))
number_subcats = number_subcats[number_subcats$step2 != "None",]

nonpancreas_order = number_subcats$step2[order(number_subcats$count[number_subcats$step2 %in% c("AdpS", "AdpV", "AdpV+S", "Liver", "Musk")])]
tissue_types = c("Adipose Only\n(Subcutaneous)", "Adipose Only\n(Visceral)", "Adipose Only\n(Both)", "Liver Only", "Muscle Only\n(Skeletal)", "Shared\nPancreas Absent", "Shared\nPancreas Present", "Pancreas Only")
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
	theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
		axis.text = element_text(size = 12))

##########################################################
# Number of colocs per trait, with and without pancreas
##########################################################

coloc_res_copy = coloc_res

loci_by_traits = coloc_res_copy %>% group_by(locus) %>% summarize(HDL = sum(grepl("High-Density-Lipoprotein", strong_coloc_traits)) > 0, WHR = sum(grepl("Waist-Hip-Ratio-BMI-Adjusted", strong_coloc_traits)) > 0, TG = sum(grepl("Triglycerides", strong_coloc_traits)) > 0, T2D = sum(grepl("Type-2-Diabetes", strong_coloc_traits)) > 0, BMI = sum(grepl("BMI", strong_coloc_traits)) > 0, FastGlu = sum(grepl("Fasting-Glucose", strong_coloc_traits)) > 0, FastInsu = sum(grepl("Fasting-Insulin", strong_coloc_traits)) > 0, ISI = sum(grepl("Insulin-Sensitivity-Index-Model-2", strong_coloc_traits)) > 0, MI = sum(grepl("M-I", strong_coloc_traits)) > 0)

loci_by_traits_with_pancreas_coloc = coloc_res_copy %>% filter(grepl("Pancreas", eqtl_file) & clpp_mod > 0.35) %>% group_by(locus) %>% summarize(HDL = sum(grepl("High-Density-Lipoprotein", strong_coloc_traits)) > 0, WHR = sum(grepl("Waist-Hip-Ratio-BMI-Adjusted", strong_coloc_traits)) > 0, TG = sum(grepl("Triglycerides", strong_coloc_traits)) > 0, T2D = sum(grepl("Type-2-Diabetes", strong_coloc_traits)) > 0, BMI = sum(grepl("BMI", strong_coloc_traits)) > 0, FastGlu = sum(grepl("Fasting-Glucose", strong_coloc_traits)) > 0, FastInsu = sum(grepl("Fasting-Insulin", strong_coloc_traits)) > 0, ISI = sum(grepl("Insulin-Sensitivity-Index-Model-2", strong_coloc_traits)) > 0, MI = sum(grepl("M-I", strong_coloc_traits)) > 0)

total_coloc_loci_per_trait = colSums(loci_by_traits[,-1])
pancreas_coloc_loci_per_trait = colSums(loci_by_traits_with_pancreas_coloc[,-1])
non_pancreas_coloc_loci_per_trait = total_coloc_loci_per_trait - pancreas_coloc_loci_per_trait

data = data.frame(list(traits = rep(names(total_coloc_loci_per_trait), 2), pancreas = c(rep("absent", 9), rep("present", 9)), loci=c(non_pancreas_coloc_loci_per_trait, pancreas_coloc_loci_per_trait)))

data$traits = factor(data$traits, levels=rev(c("ISI", "MI", "FastGlu", "FastInsu", "T2D", "WHR", "HDL", "TG", "BMI")), labels=rev(c("ISI", "M/I", "FastGlu", "FastInsu", "T2D", "WHR", "HDL", "TG", "BMI")))

g_c <- ggplot(data=data, aes(x=traits, y=loci, fill=pancreas)) +
	geom_bar(stat="identity", color="black", position=position_dodge())+
	geom_text(aes(label = loci, y = loci + 10, x = as.numeric(traits) + ifelse(pancreas=="present", 0.25, -0.25)))+
	scale_fill_manual(values=c("white", "grey"),
			  name = "pancreas\ncolocalization",
			  guide = guide_legend(reverse = TRUE)) +
	coord_flip() +
	ylab("# loci with colocalizations") +
	theme_minimal() +
	theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
		axis.text = element_text(size = 12))

##########################################################
# Now put them all together
##########################################################

g = plot_grid(NULL, g_b, g_c, labels = c('A', 'B', 'C'), nrow=1)

ggsave("output/post_coloc/plots/figure2/figure2.pdf", height=4, width=15)
