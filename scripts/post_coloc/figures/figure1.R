require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

source("scripts/post_coloc/figures/color_scheme.R")

# Plotting key info 

coloc_file = 'output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt'

# Flowchart

# Separation of the 

coloc_res = read.table(coloc_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
coloc_res$base_gwas_file = sapply(coloc_res$base_gwas_file, function(x) {s = unlist(as.vector(strsplit(x, "/"))); s[length(s)]})

###############################################
# Panel: flowchart
###############################################

# https://datascienceplus.com/how-to-build-a-simple-flowchart-with-r-diagrammer-package/

####################################################################################
# Panel: % of all loci falling out at various levels
# maybe 3 levels deep? All loci, all GWAS loci, all overlap loci, all overlapping?
####################################################################################

###############################################
# Panel: colocs per GWAS (prob add some #s)
###############################################

# How many colocs are there per GWAS...
# And then how many fall into various overlap categories

# First -- how many coloc loci are there for each GWAS; second, for each of these
# loci, what category does it fall in

coloc_res$gwas_trait=factor(x = coloc_res$gwas_trait, 
	   levels = c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz", "FastInsu_adjBMI_MAGIC_Europeans.txt.gz", "FastGlu_MAGIC_Europeans.txt.gz", "T2D_Mahajan_Europeans.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Spracklen_2020.txt.gz", "WHR-adj-BMI_GIANT_2018.txt.gz", "HDL_GLGC_Expanded.txt.gz", "TG_GLGC_Expanded.txt.gz", "BMI_GIANT_2018.txt.gz"),
	   labels = c("ISI_Genesis", "ISI_MAGIC", "FastInsu", "FastGlu", "T2D-Mahajan","T2D-Suzuki", "T2D-Xue", "T2D-Spracklen","WHR","HDL","TG","BMI"))

coloc_props = coloc_res %>% group_by(gwas_trait, locus) %>% 
	summarize(eqtl_only = (sum(nchar(strong_coloc_eqtl_genes) > 0) > 0) && (sum(nchar(strong_coloc_sqtl_genes) > 0) == 0),
		sqtl_only = (sum(nchar(strong_coloc_eqtl_genes) > 0) == 0) && (sum(nchar(strong_coloc_sqtl_genes) > 0) > 0),
		eqtl_and_sqtl = (sum(nchar(strong_coloc_eqtl_genes) > 0) > 0) && (sum(nchar(strong_coloc_sqtl_genes) > 0) > 0),
		no_colocs = (sum(nchar(strong_coloc_eqtl_genes) > 0) == 0) && (sum(nchar(strong_coloc_sqtl_genes) > 0) == 0))

# Make sure every GWAS-locus combo has exactly one category assignment 
stopifnot(sum(rowSums(coloc_props[,3:6]) != 1) == 0)

# Note: a possible problem right now is that this is out of the
# number of overlaps (coloc candidates) rather than out of the number of GWAS hits

# Might be worth actually showing both instead



gwas_qtl_types = coloc_props %>% group_by(gwas_trait) %>%
	summarize(eqtl_only = sum(eqtl_only),
		  sqtl_only = sum(sqtl_only),
		  eqtl_and_sqtl = sum(eqtl_and_sqtl),
		  no_colocs = sum(no_colocs))


mlt = melt(gwas_qtl_types)
mlt$variable = factor(mlt$variable, levels=rev(c("eqtl_only", "eqtl_and_sqtl", "sqtl_only", "no_colocs")))
colors = color_scheme[c(1,3,7,5)]

explained_ranking = rev(gwas_qtl_types$gwas_trait[order(gwas_qtl_types$no_colocs / rowSums(gwas_qtl_types[,2:5]))])
mlt$gwas_trait = factor(mlt$gwas_trait, levels = explained_ranking)

g_c = ggplot(mlt, aes(fill=variable, y=value, x=gwas_trait)) +
        geom_bar(position="fill", stat="identity", color="black") +
	scale_fill_manual(values=colors) +
	coord_flip()

# TODO: Make another version where it's out of the total number of GWAS hits instead of total number of overlapping ones

########################################################
# Panel: 
# Combined frequency histograms w/ number of candidates
# per locus or locus + number of (strong) colocs
#
# (could consider doing it also on a locus-GWAS level
# and not just locus)
########################################################

# TODO: Maybe just add in the zero bar for candidate genes as well -- some loci just don't have
# any candidate overlaps

# TODO: Maybe should really be actual count instead of density (in fact it should be, since
# total # of loci still remains the same anyway)
# but in any case, I'll move on for now; consider this a placeholder since it is not good at the moment

# A function to count the number of unique genes in a character vector, separated by semicolons
count_genes = function(x)
{
	genes = unique(strsplit(x, ";")[[1]])
	genes = genes[genes != ""]
	return(length(genes))
}

loc_sums = coloc_res %>% group_by(locus) %>% summarize(genes = sapply(paste(c(candidate_eqtl_genes, candidate_sqtl_genes), collapse=";"), count_genes), strong_colocs = sapply(paste(c(strong_coloc_eqtl_genes, strong_coloc_sqtl_genes), collapse=";"), count_genes))

points=0:100
brk = c(points-0.5, max(points)+0.5)
genes_hist = hist(loc_sums$genes, breaks=brk)
colocs_hist = hist(loc_sums$strong_colocs, breaks=brk)

plot_data = data.frame(list(count=rep(points, 2), density = c(genes_hist$density, -colocs_hist$density), type=c(rep("candidate_genes", length(points)), rep("coloc_genes", length(points)))))

plot_limit = 15
plot_data = plot_data[plot_data$count <= plot_limit,] 

g_d = ggplot(data=plot_data, aes(x=count, y=density, fill=type)) +
	geom_bar(stat = "identity", color="black") +
	scale_fill_manual(values=c("white","gray")) +
	ylim(-0.5, 0.5)


########################################################
# Now put them all together
########################################################

g = 
