require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

source("scripts/post_coloc/figures/color_scheme.R")

# Plotting key info 

coloc_file = 'output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt'
snp_gene_pairs = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_snp_gene_pairs_prefiltering.txt"

clpp_mod_threshold = 0.4

coloc_res = read.table(coloc_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
coloc_res$base_gwas_file = sapply(coloc_res$base_gwas_file, function(x) {s = unlist(as.vector(strsplit(x, "/"))); s[length(s)]})

coloc_res$tissue = coloc_res$eqtl_file
coloc_res$tissue = gsub("data/eqtls/gtex_v8/", "", coloc_res$tissue)
coloc_res$tissue = gsub("data/sqtls/gtex_v8/", "", coloc_res$tissue)
coloc_res$tissue = gsub(".sQTLs.txt.gz", "", coloc_res$tissue)
coloc_res$tissue = gsub(".allpairs.txt.gz.eQTLs.txt.gz", "", coloc_res$tissue)

overlap_table = read.table(snp_gene_pairs, header=TRUE, sep="\t", stringsAsFactors=FALSE)
overlap_table$tissue = overlap_table$lookup_file

overlap_table$tissue = gsub("data/eqtls/gtex_v8/", "", overlap_table$tissue)
overlap_table$tissue = gsub("data/sqtls/gtex_v8/", "", overlap_table$tissue)
overlap_table$tissue = gsub(".sQTLs.txt.gz", "", overlap_table$tissue)
overlap_table$tissue = gsub(".allpairs.txt.gz.eQTLs.txt.gz", "", overlap_table$tissue)

overlap_counts = overlap_table %>% group_by(tissue) %>% summarize(gene_count=length(unique(ensembl)))
test_counts = coloc_res %>% group_by(tissue) %>% summarize(gene_count=length(unique(ensembl)))

coloc_passing =  coloc_res %>% filter(clpp_mod > clpp_mod_threshold)
strong_coloc_counts = coloc_passing %>% group_by(tissue) %>% summarize(gene_count=length(unique(ensembl)))


###############################################
# Panel: flowchart
###############################################

# Dummy data; load real data later
tissue_colors = c("#FF6600", "#FFAA00", "#AABB66", "#995522", "#AAAAFF")
tissues = c("Adipose-Subcutaneous", "Adipose-Visceral", "Liver", "Pancreas", "Skeletal-Muscle")
#gene_overlaps = c(50000, 50000, 50000, 50000, 50000)
gene_overlaps = overlap_counts$gene_count
#gene_qtl_overlaps = c(3000, 6000, 2000, 1000, 2000)
gene_qtl_overlaps = test_counts$gene_count
#gene_colocs = c(100, 200, 300, 500, 200)
gene_colocs = strong_coloc_counts$gene_count

# Do a bit of math to figure out the geometry of the plot

# Specify relative heights of the different elements to be drawn in the picture
bar_heights=c(1,1,1)
connector_heights=c(0.7,0.7)
bar_spacing=c(0.4,0.1,0.1,0.1,0.1,0.4)

# Might not be necessary
element_heights = c(bar_spacing[1], bar_heights[1], bar_spacing[2], connector_heights[1], bar_spacing[3], bar_heights[2], bar_spacing[4], connector_heights[2], bar_spacing[5], bar_heights[3], bar_spacing[6])
start_heights = cumsum(c(0,element_heights))

#overlap_percents = gene_overlaps / sum(gene_overlaps)
#total_qtl_overlap_percents = gene_qtl_overlaps / sum(gene_overlaps)
#total_coloc_percents = gene_colocs / sum(gene_overlaps)

#qtl_overlap_percents = gene_qtl_overlaps / sum(gene_qtl_overlaps)
#coloc_percents = gene_colocs / sum(gene_colocs)


# Assemble coordinates of gene overlap boxes
gene_overlap_data = draw_gene_overlap_boxes(start_heights, gene_overlaps, tissues)
gene_qtl_overlap_data = draw_gene_qtl_overlap_boxes(start_heights, gene_qtl_overlaps, tissues)
gene_coloc_data = draw_gene_coloc_boxes(start_heights, gene_colocs, tissues)

first_connector_data = draw_first_connectors(start_heights, gene_overlaps, gene_qtl_overlaps, tissues)
second_connector_data = draw_second_connectors(start_heights, gene_qtl_overlaps, gene_colocs, tissues)


unique_genes = length(unique(overlap_table$ensembl))
unique_loci = length(unique(overlap_table$locus))

unique_qtl_genes = length(unique(coloc_res$ensembl))
unique_qtl_loci = length(unique(coloc_res$locus))

unique_coloc_genes = length(unique(coloc_passing$ensembl))
unique_coloc_loci = length(unique(coloc_passing$locus))

label_x = -sum(gene_overlaps) * 3 / 4

g_a = ggplot() + theme_void() +
	geom_polygon(data=gene_overlap_data, mapping=aes(x=x, y=y, group=groups, fill=groups)) + 
	geom_polygon(data=gene_qtl_overlap_data, mapping=aes(x=x, y=y, group=groups, fill=groups)) + 
	geom_polygon(data=gene_coloc_data, mapping=aes(x=x, y=y, group=groups, fill=groups)) + 
	geom_polygon(data=first_connector_data, mapping=aes(x=x, y=y, group=groups, fill=groups), alpha=0.5) + 
	geom_polygon(data=second_connector_data, mapping=aes(x=x, y=y, group=groups, fill=groups), alpha=0.5) + 
	geom_text(aes(label_x, -(start_heights[2] + start_heights[3]) / 2, label=sprintf("gene overlaps\n\n%d unique genes\n%d unique loci", unique_genes, unique_loci)), hjust=0) + 
	geom_text(aes(label_x, -(start_heights[6] + start_heights[7]) / 2, label=sprintf("significant QTL overlaps\n\n%d unique genes\n%d unique loci", unique_qtl_genes, unique_qtl_loci)), hjust=0) + 
	geom_text(aes(label_x, -(start_heights[10] + start_heights[11]) / 2, label=sprintf("colocalizations\n\n%d unique genes\n%d unique loci", unique_coloc_genes, unique_coloc_loci)), hjust=0) +
	scale_fill_manual(values=tissue_colors) +
	theme(legend.title=element_blank())
	#geom_text(aes(1,0,label='N/A')) + 
	#geom_text(aes(1,1,label='N/A')) + 
	#geom_text(aes(0,1,label='N/A'))

# To add as labels: # unique genes at the level, # unique loci


g_a
# https://datascienceplus.com/how-to-build-a-simple-flowchart-with-r-diagrammer-package/

####################################################################################
# Panel: % of all loci falling out at various levels
# maybe 3 levels deep? All loci, all GWAS loci, all overlap loci, all overlapping?
####################################################################################

g_b = NULL

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
mlt$variable = factor(mlt$variable, labels=rev(c("eqtl only", "eqtl and sqtl", "sqtl only", "no colocs")), levels=rev(c("eqtl_only", "eqtl_and_sqtl", "sqtl_only", "no_colocs")))
colors = color_scheme[c(1,3,7,5)]

explained_ranking = rev(gwas_qtl_types$gwas_trait[order(gwas_qtl_types$no_colocs / rowSums(gwas_qtl_types[,2:5]))])
mlt$gwas_trait = factor(mlt$gwas_trait, levels = explained_ranking)

g_c = ggplot(mlt, aes(fill=variable, y=value, x=gwas_trait)) +
        geom_bar(position="fill", stat="identity", color="black") +
	scale_fill_manual(values=colors) +
	theme_minimal() +
	theme(axis.title.y = element_blank(), legend.title=element_blank()) +
	coord_flip() +
	ylab("% of tested loci") 

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

plot_data = data.frame(list(count=rep(points, 2), density = c(genes_hist$density, -colocs_hist$density), type=c(rep("candidate genes", length(points)), rep("colocalized genes", length(points)))))

plot_limit = 15
plot_data = plot_data[plot_data$count <= plot_limit,] 

g_d = ggplot(data=plot_data, aes(x=count, y=density, fill=type)) +
	geom_bar(stat = "identity", color="black") +
	scale_fill_manual(values=c("white","gray")) +
	ylim(-0.5, 0.5) +
	theme_minimal() +
	theme(axis.title.y = element_blank(), legend.title=element_blank())


########################################################
# Now put them all together
########################################################

first_row = plot_grid(g_a, labels = c('A'), nrow=1)
second_row = plot_grid(g_b, g_c, g_d, labels = c('B', 'C', 'D'), nrow = 3)

gg_all = plot_grid(first_row, second_row, labels=c('', ''), ncol=2)


##########################################
# Auxiliary functions
##########################################

draw_gene_overlap_boxes = function(start_heights, gene_overlaps, tissues)
{
	vertices_per_shape = 4
	shape_count = length(gene_overlaps)
	vertex_count = shape_count * vertices_per_shape

	x = rep(0, vertex_count)
	y = rep(0, vertex_count)
	groups = rep("", vertex_count)

	x_scaffold = cumsum(c(0, gene_overlaps))
	x_scaffold = x_scaffold - mean(x_scaffold)

	for (i in 1:shape_count)
	{
		x1 = x_scaffold[i]
		x2 = x_scaffold[i+1]
		y1 = start_heights[2] 
		y2 = start_heights[3] 

		groups[((i-1)*vertices_per_shape+1):(i*vertices_per_shape)] = as.character(tissues[i])
		x[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+4)] = x1
		x[c((i-1)*(vertices_per_shape)+2, (i-1)*(vertices_per_shape)+3)] = x2	
		y[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+2)] = y1	
		y[c((i-1)*(vertices_per_shape)+3, (i-1)*(vertices_per_shape)+4)] = y2	
	}
	return(data.frame(list(x=x, y=-y, groups=groups)))
}

draw_gene_qtl_overlap_boxes = function(start_heights, gene_qtl_overlaps, tissues)
{
	vertices_per_shape = 4
	shape_count = length(gene_overlaps)
	vertex_count = shape_count * vertices_per_shape

	x = rep(0, vertex_count)
	y = rep(0, vertex_count)
	groups = rep("", vertex_count)

	x_scaffold = cumsum(c(0, gene_qtl_overlaps))
	x_scaffold = x_scaffold - mean(x_scaffold)

	for (i in 1:shape_count)
	{
		x1 = x_scaffold[i]
		x2 = x_scaffold[i+1]
		y1 = start_heights[6] 
		y2 = start_heights[7] 

		groups[((i-1)*vertices_per_shape+1):(i*vertices_per_shape)] = as.character(tissues[i])
		x[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+4)] = x1
		x[c((i-1)*(vertices_per_shape)+2, (i-1)*(vertices_per_shape)+3)] = x2	
		y[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+2)] = y1	
		y[c((i-1)*(vertices_per_shape)+3, (i-1)*(vertices_per_shape)+4)] = y2	
	}
	return(data.frame(list(x=x, y=-y, groups=groups)))
}

draw_gene_coloc_boxes = function(start_heights, gene_colocs, tissues)
{
	# Assemble coordinates of coloc count boxes
	vertices_per_shape = 4
	shape_count = length(gene_overlaps)
	vertex_count = shape_count * vertices_per_shape

	x = rep(0, vertex_count)
	y = rep(0, vertex_count)
	groups = rep("", vertex_count)

	x_scaffold = cumsum(c(0, gene_colocs))
	x_scaffold = x_scaffold - mean(x_scaffold)

	for (i in 1:shape_count)
	{
		x1 = x_scaffold[i]
		x2 = x_scaffold[i+1]
		y1 = start_heights[10] 
		y2 = start_heights[11] 

		groups[((i-1)*vertices_per_shape+1):(i*vertices_per_shape)] = as.character(tissues[i])
		x[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+4)] = x1
		x[c((i-1)*(vertices_per_shape)+2, (i-1)*(vertices_per_shape)+3)] = x2	
		y[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+2)] = y1	
		y[c((i-1)*(vertices_per_shape)+3, (i-1)*(vertices_per_shape)+4)] = y2	
	}
	return(data.frame(list(x=x, y=-y, groups=groups)))
}

draw_first_connectors = function(start_heights, gene_overlaps, gene_qtl_overlaps, tissues)
{
	# Assemble coordinates of coloc count boxes
	vertices_per_shape = 4
	shape_count = length(gene_overlaps)
	vertex_count = shape_count * vertices_per_shape

	x = rep(0, vertex_count)
	y = rep(0, vertex_count)
	groups = rep("", vertex_count)

	x_scaffold = cumsum(c(0, gene_overlaps))
	x_scaffold = x_scaffold - mean(x_scaffold)
	x_scaffold2 = cumsum(c(0, gene_qtl_overlaps))
	x_scaffold2 = x_scaffold2 - mean(x_scaffold2)

	for (i in 1:shape_count)
	{
		x1 = x_scaffold[i]
		x2 = x_scaffold[i+1]
		x3 = x_scaffold2[i+1]
		x4 = x_scaffold2[i]
		y1 = start_heights[4] 
		y2 = start_heights[5] 

		groups[((i-1)*vertices_per_shape+1):(i*vertices_per_shape)] = as.character(tissues[i])
		x[c((i-1)*(vertices_per_shape)+1)] = x1
		x[c((i-1)*(vertices_per_shape)+2)] = x2	
		x[c((i-1)*(vertices_per_shape)+3)] = x3
		x[c((i-1)*(vertices_per_shape)+4)] = x4	
		y[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+2)] = y1	
		y[c((i-1)*(vertices_per_shape)+3, (i-1)*(vertices_per_shape)+4)] = y2	
	}
	return(data.frame(list(x=x, y=-y, groups=groups)))

}

draw_second_connectors = function(start_heights, gene_qtl_overlaps, gene_colocs, tissues)
{
	# Assemble coordinates of coloc count boxes
	vertices_per_shape = 4
	shape_count = length(gene_overlaps)
	vertex_count = shape_count * vertices_per_shape

	x = rep(0, vertex_count)
	y = rep(0, vertex_count)
	groups = rep("", vertex_count)

	x_scaffold = cumsum(c(0, gene_qtl_overlaps))
	x_scaffold = x_scaffold - mean(x_scaffold)
	x_scaffold2 = cumsum(c(0, gene_colocs))
	x_scaffold2 = x_scaffold2 - mean(x_scaffold2)

	for (i in 1:shape_count)
	{
		x1 = x_scaffold[i]
		x2 = x_scaffold[i+1]
		x3 = x_scaffold2[i+1]
		x4 = x_scaffold2[i]
		y1 = start_heights[8] 
		y2 = start_heights[9] 

		groups[((i-1)*vertices_per_shape+1):(i*vertices_per_shape)] = as.character(tissues[i])
		x[c((i-1)*(vertices_per_shape)+1)] = x1
		x[c((i-1)*(vertices_per_shape)+2)] = x2	
		x[c((i-1)*(vertices_per_shape)+3)] = x3
		x[c((i-1)*(vertices_per_shape)+4)] = x4	
		y[c((i-1)*(vertices_per_shape)+1, (i-1)*(vertices_per_shape)+2)] = y1	
		y[c((i-1)*(vertices_per_shape)+3, (i-1)*(vertices_per_shape)+4)] = y2	
	}
	return(data.frame(list(x=x, y=-y, groups=groups)))


}
