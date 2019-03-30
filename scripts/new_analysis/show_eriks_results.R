# Author: Mike Gloudemans

require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)

# Define the list of traits that we actually care about.
# We will load all traits but then subset down to these.
kept_traits = c("/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_FastInsu_adjBMI_MAGIC_Europeans.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_FastGlu_MAGIC_Europeans.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_BMI_GIANT_Europeans.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_WHRadjBMI_GIANT_Europeans.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_T2D_DIAGRAM_European.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_TG_GCLC_Mixed.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_HDL_GCLC_Mixed.txt.gz",
		"/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_MI_adjBMI.txt.gz")

#
# Load and prep COLOC results files
#

# Load all files into R for analysis
folders = dir("/users/mgloud/projects/brain_gwas/output/insulin-full-analysis/")
tabs = list()
i = 1
for (folder in folders)
{
	files = dir(paste0("/users/mgloud/projects/brain_gwas/output/insulin-full-analysis/", folder))
	files = files[grep("clpp" ,files)]
	for (file in files)
	{
		tabs[[i]] = read.table(paste("/users/mgloud/projects/brain_gwas/output/insulin-full-analysis", folder, file, sep="/"), header=TRUE)
		i = i + 1
	}
}
results = do.call(rbind, tabs)

# A bit of stuff just for more understandable labels
results$ensembl = sapply(as.character(results$feature), function(x) {strsplit(x, "\\.")[[1]][1]})
results$gwas_trait = gsub(".prepared.txt.gz", "", results$gwas_trait)
results$gwas_trait = gsub("_Mixed", "", results$gwas_trait)
results$gwas_trait = gsub("_Europeans", "", results$gwas_trait)
results$gwas_trait = gsub("_European", "", results$gwas_trait)
results$gwas_trait = gsub("_AllSNPs", "", results$gwas_trait)
results$gwas_trait = gsub("_MAGIC", "", results$gwas_trait)
results$gwas_trait = gsub("_GCLC", "", results$gwas_trait)
results$gwas_trait = gsub("_DIAGRAM", "", results$gwas_trait)
results$gwas_trait = gsub("_GIANT", "", results$gwas_trait)
results$gwas_trait = gsub("_CARDIoGRAMplusC4D", "", results$gwas_trait)
results$gwas_trait = gsub(".txt", "", results$gwas_trait)
results$gwas_trait = gsub("/users/mgloud/projects/gwas/data/prepared/", "", results$gwas_trait)
results$gwas_trait = gsub("/users/mgloud/projects/gwas/data/GENESIS/formatted/", "", results$gwas_trait)
results$gwas_trait = gsub("GWAS_", "", results$gwas_trait)
results$gwas_trait = gsub(".txt.gz", "", results$gwas_trait)
results$eqtl_file = gsub("_allpairs_txt_gz", "", results$eqtl_file)

# Get HGNC gene names
genes = read.table("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/mart_export.txt", header=TRUE, sep="\t")
genes = genes[,-2]
names(genes) = c("ensembl", "hgnc")
genes = genes[!duplicated(genes$ensembl),]
results = merge(genes, results, all.y=TRUE)

# Output a file in Bruna's requested format

# Load rsIDs
# TODO: 
# finemap = merge(rsids, finemap)

ensembl_locs = read_delim("/srv/scratch/restricted/rare_diseases/data/annot/gencode.v19.annotation.genes.bed.gz", delim="\t", col_names=FALSE)
colnames(ensembl_locs) = c("chrom", "gene_start", "gene_end", "feature")

results$gene_start = ensembl_locs[match(results$feature, ensembl_locs$feature),]$gene_start
results$gene_end = ensembl_locs[match(results$feature, ensembl_locs$feature),]$gene_end

results$chr = sapply(as.character(results$ref_snp), function(x) {strsplit(x, "_")[[1]][1]})
results$pos = sapply(as.character(results$ref_snp), function(x) {strsplit(x, "_")[[1]][2]})

results$min_gwas_pval = 10^-(results$X.log_gwas_pval)
results$min_eqtl_pval = 10^-(results$X.log_eqtl_pval)

# Cluster by locus
ids = unique(as.character(results$ref_snp))
ids = ids[order(ids)]
chr = sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][1]})
pos = as.numeric(sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][2]}))
loc_nums = rep(0, length(ids))
loc_nums[1] = 1
for (i in 2:length(ids))
{
	# Check if there's a SNP above in the list within 1 MB of this SNP]
	same = (chr[1:(i-1)] == chr[i]) & (abs(pos[1:(i-1)] - pos[i]) < 1000000)
	if (length(which(same)) != 0)
	{
		loc_nums[i] = loc_nums[which(same)[1]]
	}
	else
	{
		loc_nums[i] = max(loc_nums) + 1
	}
}

results$loc_num = sapply(results$ref_snp, function(x) {
		loc_nums[which(ids == x)]
	})
# Just check the order to make sure it looks good
sub = results[!duplicated(results$ref_snp),]
sub[order(as.character(sub$ref_snp)),][c("ref_snp", "loc_num")]

results_output = results[c("chr", "pos", "loc_num", "min_gwas_pval", "min_eqtl_pval", "gwas_trait", "eqtl_file", "ensembl", "hgnc", "gene_start", "gene_end", "clpp", "clpp_mod", "n_snps")]
write.table(results_output, "/users/mgloud/projects/insulin_resistance/output/clpp_results.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


# Initalize a 3D CLPP array of tissue x gene x gwas_trait
# Maybe initialize another indicator array to show which ones meet all the cutoffs (could then remove slices with no results meeting the cutoffs)
tissues = unique(results$eqtl_file)
genes = unique(results$feature)
gwas_traits = unique(results$gwas_trait)

#min_clpp = 0.001
min_clpp = 0
clpp = array(min_clpp, dim = c(length(tissues), length(genes), length(gwas_traits)))

dimnames(clpp)[[1]] = tissues
dimnames(clpp)[[2]] = genes
dimnames(clpp)[[3]] = gwas_traits

# Fill the array
for (i in 1:dim(clpp)[1])
{
	for (j in 1:dim(clpp)[2])
	{
		for (k in 1:dim(clpp)[3])
		{
			clpp[i,j,k] = max(results[(results$eqtl_file == tissues[i]) & (results$feature == genes[j]) & (results$gwas_trait == gwas_traits[k]),]$clpp_mod, clpp[i,j,k])
		}
	}
}

clpp = clpp[,apply(clpp, 2, max) >= 0.3,]


# Cluster the array within each plane
# May need to read a bit about how to do this and how to choose correlation parameters...but don't waste too much time on it
# TODO: (I'll skip it for now, just focus on consistency)

# Plot it in interesting ways? For now though just do layered heatmaps, one for each tissue
# (Actually, would be nice to be able to view all 3 possible slice orientations)

# Make one plot for each tissue
for (i in 1:dim(clpp)[1])
{
	tissue = tissues[i]

	my_matrix = as.data.frame(as.matrix(clpp[i,,]))
	my_matrix$gene = rownames(my_matrix)

	# Plot heatmap
	heat = melt(my_matrix)
	heat$plot_vals = heat$value
	heat$binned_vals = NA
	heat$binned_vals[heat$plot_vals <= 1] = 8
	heat$binned_vals[heat$plot_vals < 0.9] = 7
	heat$binned_vals[heat$plot_vals < 0.8] = 6
	heat$binned_vals[heat$plot_vals < 0.7] = 5
	heat$binned_vals[heat$plot_vals < 0.6] = 4
	heat$binned_vals[heat$plot_vals < 0.5] = 3
	heat$binned_vals[heat$plot_vals < 0.4] = 2
	heat$binned_vals[heat$plot_vals < 0.3] = 1

	heat$binned_vals = as.factor(heat$binned_vals)
	
	heat$gene = factor(heat$gene, levels=unique(heat$gene))
	heat$variable = factor(heat$variable, levels=unique(heat$variable))

	heat$strike = ""
	heat$strike[is.na(heat$binned_vals)] = "X"
	heat$binned_vals[is.na(heat$binned_vals)] = 1

	g = ggplot(data = heat, aes(x = variable, y = gene)) +
			geom_tile(aes(fill = binned_vals), color="grey80", size=0.5) +
			scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(8), name="CLPP_mod", labels = c("0-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"),na.value = "white") +
			theme(axis.text.y=element_text(size=10),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
			labs(y = "Nearest Gene --> Colocalized Gene") +
			labs(x = "GWAS Trait") + 
			geom_text(label=heat$strike, size=4, color="red", alpha=0.4)

	ggsave(paste0("/users/mgloud/projects/insulin_resistance/output/new/insulin_", tissue, "_full_clpp.png"), width = 8, height = 30, units = "in", dpi = 300, limitsize=FALSE)
}


# Make one plot for each gene
for (i in 1:dim(clpp)[2])
{
	gene = genes[i]

	my_matrix = as.data.frame(as.matrix(clpp[,i,]))
	my_matrix$tissue = rownames(my_matrix)

	# Plot heatmap
	heat = melt(my_matrix)
	heat$plot_vals = heat$value
	heat$binned_vals = NA
	heat$binned_vals[heat$plot_vals <= 1] = 8
	heat$binned_vals[heat$plot_vals < 0.9] = 7
	heat$binned_vals[heat$plot_vals < 0.8] = 6
	heat$binned_vals[heat$plot_vals < 0.7] = 5
	heat$binned_vals[heat$plot_vals < 0.6] = 4
	heat$binned_vals[heat$plot_vals < 0.5] = 3
	heat$binned_vals[heat$plot_vals < 0.4] = 2
	heat$binned_vals[heat$plot_vals < 0.3] = 1

	heat$binned_vals = as.factor(heat$binned_vals)
	
	heat$tissue = factor(heat$tissue, levels=unique(heat$tissue))
	heat$variable = factor(heat$variable, levels=unique(heat$variable))

	heat$strike = ""
	heat$strike[is.na(heat$binned_vals)] = "X"
	heat$binned_vals[is.na(heat$binned_vals)] = 1

	g = ggplot(data = heat, aes(x = variable, y = tissue)) +
			geom_tile(aes(fill = binned_vals), color="grey80", size=0.5) +
			scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(8), name="CLPP_mod", labels = c("0-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"),na.value = "white") +
			theme(axis.text.y=element_text(size=10),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
			labs(y = "Tissue") +
			labs(x = "GWAS Trait") + 
			geom_text(label=heat$strike, size=4, color="red", alpha=0.4)

	ggsave(paste0("/users/mgloud/projects/insulin_resistance/output/new/insulin_", gene, "_full_clpp.png"), width = 8, height = 30, units = "in", dpi = 300, limitsize=FALSE)
}


# Make one plot for each GWAS
for (i in 1:dim(clpp)[3])
{
	gwas= gwas_traits[i]

	my_matrix = as.data.frame(t(as.matrix(clpp[,,i])))
	my_matrix$gene = rownames(my_matrix)

	# Plot heatmap
	heat = melt(my_matrix)
	heat$plot_vals = heat$value
	heat$binned_vals = NA
	heat$binned_vals[heat$plot_vals <= 1] = 8
	heat$binned_vals[heat$plot_vals < 0.9] = 7
	heat$binned_vals[heat$plot_vals < 0.8] = 6
	heat$binned_vals[heat$plot_vals < 0.7] = 5
	heat$binned_vals[heat$plot_vals < 0.6] = 4
	heat$binned_vals[heat$plot_vals < 0.5] = 3
	heat$binned_vals[heat$plot_vals < 0.4] = 2
	heat$binned_vals[heat$plot_vals < 0.3] = 1

	heat$binned_vals = as.factor(heat$binned_vals)
	
	heat$gene = factor(heat$gene, levels=unique(heat$gene))
	heat$variable = factor(heat$variable, levels=unique(heat$variable))

	heat$strike = ""
	heat$strike[is.na(heat$binned_vals)] = "X"
	heat$binned_vals[is.na(heat$binned_vals)] = 1

	g = ggplot(data = heat, aes(x = variable, y = gene)) +
			geom_tile(aes(fill = binned_vals), color="grey80", size=0.5) +
			scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(8), name="CLPP_mod", labels = c("0-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"),na.value = "white") +
			theme(axis.text.y=element_text(size=10),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
			labs(y = "Gene") +
			labs(x = "Tissue") + 
			geom_text(label=heat$strike, size=4, color="red", alpha=0.4)

	ggsave(paste0("/users/mgloud/projects/insulin_resistance/output/new/insulin_", gwas, "_full_clpp.png"), width = 8, height = 30, units = "in", dpi = 300, limitsize=FALSE)
}


stopifnot(FALSE)


# Use 3D heatmap vis? We'll see, this isn't top priority right now.
# A good long-term goal
