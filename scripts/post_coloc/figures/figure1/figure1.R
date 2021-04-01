require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

#################################
# Constants and settings
#################################

source("scripts/post_coloc/figures/color_scheme.R")

coloc_file = 'output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt'
snp_gene_pairs = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_snp_gene_pairs_prefiltering.txt"
gwas_snps_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/all_gwas_snps.txt"
all_loci_file = "data/ldetect/fourier_ls-all.hg38.connected.bed"


clpp_mod_threshold = 0.4

main = function() {

	#################################
	# Bit of data preprocessing
	#################################

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
	sqtl_coloc_passing =  coloc_passing %>% filter(grepl("sQTL", eqtl_short))
	sqtl_strong_coloc_counts = sqtl_coloc_passing %>% group_by(tissue) %>% summarize(gene_count=length(unique(ensembl)))
	eqtl_coloc_passing =  coloc_passing %>% filter(grepl("eQTL", eqtl_short))
	eqtl_strong_coloc_counts = eqtl_coloc_passing %>% group_by(tissue) %>% summarize(gene_count=length(unique(ensembl)))

	gwas_hits = read.table(gwas_snps_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

	all_genome_loci = read.table(all_loci_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

	# Dummy data; load real data later
	tissue_colors = c("#FF6600", "#FFAA00", "#AABB66", "#AAAAFF", "#995522")
	tissues = c("Adipose (Sub)", "Adipose (Visc)", "Liver", "Muscle", "Pancreas")

	###############################################
	# Panel: flowchart
	###############################################

	# This figure was done manually in inkscape; see `figure1a-metrics` for the
	# code used to get the numbers themselves


	####################################################################################
	# Panel: comparison of GTEx tissue sample size and # colocs found in that tissue
	####################################################################################

	sample_sizes = c(  581, # Adipose_Subcutaneous 
			   469, # Adipose_Visceral_Omentum
			   208, # Liver
			   706, # Muscle_Skeletal
			   305) # Pancreas

	egenes = c(        15607, # Adipose_Subcutaneous 
			   12482, # Adipose_Visceral_Omentum
			   5734, # Liver
			   13532, # Muscle_Skeletal
			   9660) # Pancreas
	
	sgenes = c(        5113, # Adipose_Subcutaneous 
			   4210, # Adipose_Visceral_Omentum
			   1485, # Liver
			   4056, # Muscle_Skeletal
			   2250) # Pancreas


	# (these numbers were generated above in the section for the panel g_a that is no longer used)
	tissue_qtls = data.frame(list(tissue=strong_coloc_counts$tissue, eqtl=eqtl_strong_coloc_counts$gene_count, sqtl=sqtl_strong_coloc_counts$gene_count, either=strong_coloc_counts$gene_count, sample_sizes = sample_sizes , egenes=egenes, sgenes=sgenes))

	tissue_eqtls = data.frame(list(x=tissue_qtls$egenes, y=tissue_qtls$eqtl, type="eqtl", tissue=tissues))
	tissue_sqtls = data.frame(list(x=tissue_qtls$sgenes, y=tissue_qtls$sqtl, type="sqtl", tissue=tissues))

	all = rbind(tissue_eqtls, tissue_sqtls)

	
	# This part is supplementary and was ultimately not included in the paper (was not ultimately included)
	g_b_supp = ggplot(all, aes(x=x, y=y, color=tissue, shape=type)) +
		geom_point(size=3) +
		theme_minimal() +
		scale_color_manual(values=tissue_colors)+
		xlim(c(0,max(all$x)))+
		ylim(c(0,max(all$y)))+
		xlab("# eGenes / sGenes")+
		ylab("# eQTL / sQTL tests colocalized")
	g_b_supp

	# This is the part for the main figure panel 1B

	tissue_eqtls = data.frame(list(x=tissue_qtls$sample_sizes, y=tissue_qtls$eqtl, type="eqtl", tissue=tissues))
	tissue_sqtls = data.frame(list(x=tissue_qtls$sample_sizes, y=tissue_qtls$sqtl, type="sqtl", tissue=tissues))
	
	all = rbind(tissue_eqtls, tissue_sqtls)

	g_b = ggplot(all, aes(x=x, y=y, color=tissue, shape=type)) +
		geom_point(size=3) +
		theme_minimal() +
		scale_color_manual(values=tissue_colors)+
		xlim(c(0,max(all$x)))+
		ylim(c(0,max(all$y)))+
		xlab("Tissue sample size")+
		ylab("# eQTL / sQTL tests colocalized")
	g_b


	####################################
	# eQTL / sQTL sharing 
	####################################

	# See separate script for counting; the figure itself was generated
	# manually in Inkscape

	########################################################
	# Panel: 
	# Combined frequency of number of candidates
	# per locus or locus + number of (strong) colocs
	########################################################

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


	loc_sums$genes = pmin(loc_sums$genes, 20)
	loc_sums$strong_colocs = pmin(loc_sums$strong_colocs, 20)

	g_d = ggplot(data=loc_sums, aes(x=genes+rnorm(length(genes),0,0.1), y=strong_colocs+rnorm(length(genes),0,0.1))) +
		theme_minimal() +
		geom_point(color="red",stroke=0,size=2,alpha=0.2)+
		xlab("Number of candidate genes")+
		ylab("Number of colocalized genes")

	########################################################
	# Now put them all together
	########################################################

	first_row = plot_grid(NULL, labels = c('A'), nrow=1)
	second_row = plot_grid(g_b, NULL, g_d, labels = c('B', 'C', 'D'), ncol=3)

	gg_all = plot_grid(first_row, second_row, labels=c('', ''), rel_heights = c(3,2), nrow=2)

	gg_all
	ggsave("output/post_coloc/plots/figure1/figure1.pdf", height=8, width=12)

}

main()
