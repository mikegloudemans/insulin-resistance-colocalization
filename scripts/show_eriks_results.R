# Author: Mike Gloudemans
# Date created: 4/10/2018

require(reshape)
require(ggplot2)
require(dplyr)

#
# Load files
#

files = c(
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/BMI_GIANT_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/CHD_CARDIoGRAMplusC4D_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/FastGlu_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/FastInsu_adjBMI_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
#	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/GWAS_GENESIS_MI_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/HDL_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/LDL_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/T2D_DIAGRAM_European_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/TG_GCLC_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
	  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-10-19_13-29-57.409767_insulin_resistance_perturbations_liver_only/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt"
	  )

# Finemap results
finemap_results = list()
for (i in 1:length(files))
{
	finemap_results[[i]] = read.table(files[i], header=TRUE)
}
finemap = do.call(rbind, finemap_results)
finemap$ensembl = sapply(as.character(finemap$feature), function(x) {strsplit(x, "\\.")[[1]][1]})
finemap$gwas_trait = gsub(".prepared.txt.gz", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_Mixed", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_Europeans", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_European", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_AllSNPs", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_MAGIC", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_GCLC", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_DIAGRAM", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_GIANT", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("_CARDIoGRAMplusC4D", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("/users/mgloud/projects/gwas/data/prepared/", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("/users/mgloud/projects/gwas/data/GENESIS/formatted/", "", finemap$gwas_trait)
finemap$gwas_trait = gsub("GWAS_", "", finemap$gwas_trait)
finemap$gwas_trait = gsub(".txt.gz", "", finemap$gwas_trait)

finemap$eqtl_file = sapply(as.character(finemap$eqtl_file), function(x)
       {
               s = strsplit(x, "_")[[1]]
               return(paste(s[1:(length(s)-3)], collapse="_"))
       })

# Load rsIDs
rsids = read.table("/users/mgloud/projects/insulin_resistance/output/filtered_insulin_snps_to_test.txt", fill=TRUE, header=TRUE)
rsids = rsids[,1:3]
rsids$ref_snp = paste(rsids[,2], rsids[,3], sep="_")
rsids = rsids[,-c(2,3)]
colnames(rsids)[1] = "rsid"

# Get HGNC gene names
genes = read.table("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/mart_export.txt", header=TRUE, sep="\t")
genes = genes[,-2]
names(genes) = c("ensembl", "hgnc")

genes = genes[!duplicated(genes$ensembl),]

finemap = merge(rsids, finemap)
finemap = merge(genes, finemap)

require(readr)
ensembl_locs = read_delim("/srv/scratch/restricted/rare_diseases/data/annot/gencode.v19.annotation.genes.bed.gz", delim="\t", col_names=FALSE)
colnames(ensembl_locs) = c("chrom", "gene_start", "gene_end", "feature")

finemap$gene_start = ensembl_locs[match(finemap$feature, ensembl_locs$feature),]$gene_start
finemap$gene_end = ensembl_locs[match(finemap$feature, ensembl_locs$feature),]$gene_end

finemap$chr = sapply(finemap$ref_snp, function(x) {strsplit(x, "_")[[1]][1]})
finemap$pos = sapply(finemap$ref_snp, function(x) {strsplit(x, "_")[[1]][2]})

results_output = finemap[c("rsid", "chr", "pos", "gwas_trait", "X.log_gwas_pval", "feature", "hgnc", "gene_start", "gene_end", "X.log_eqtl_pval", "clpp", "clpp_mod")]
write.table(results_output, "/users/mgloud/projects/insulin_resistance/output/clpp_results.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

#
# First, get the list of SNP/gene pairs passing the desired cutoff, so that we can
# show them in individual plots
#
dim(finemap)	# Total number of tests performed
pass_pval_cutoffs = finemap[finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 5,]
dim(pass_pval_cutoffs)	# Total number of tests (gwas SNP - eqtl gene pairs) passing eQTL and GWAS pval thresholds
final_set = pass_pval_cutoffs[pass_pval_cutoffs$clpp > 0.01,]
dim(final_set)	# Total number of tests passing our colocalization criteria

# Remove duplicates from this list
final_set = final_set[!duplicated(paste(final_set$ref_snp, final_set$feature, sep="_")),]

# Produce a list of genes for prioritization in follow-up studies, ranked by CLPP_mod score.
# Order by max CLPP score
# Note that clusters haven't been removed yet from this matrix
priority = final_set[rev(order(final_set$clpp)),]
priority = priority[!duplicated(priority$feature),] # Fine to do this now, since it guarantees we're taking the top coloc at least

#
print(dim(priority)[1])
# Number of SNPs that had at least one hit
print(length(unique(priority$ref_snp)))
# Contrast with total number of SNPs
print(length(unique(finemap$ref_snp)))

# Print a table showing SNP, gene, CLPP score in top context, coloc rank among all genes associated with that SNP.
# Other info can be found from the heatmaps

# We have most of this info already, but we need the column showing rank of this gene among all colocalized
# genes at this locus
priority$rank_at_locus = sapply(1:dim(priority)[1], function(x) 
	   	{
			return(sum(priority$ref_snp[1:x] == priority$ref_snp[x]))
		}
	   )

priority$ensembl = as.character(priority$ensembl)
priority$hgnc = as.character(priority$hgnc)
priority$hgnc[priority$hgnc == ""] = priority$ensembl[priority$hgnc == ""]

# Then can plot the total number of colocalizations found for each gene and include that in the presentation too
last_only = priority[!rev(duplicated(rev(priority$rsid))),]	# Get last appearance of each rsid in the list
hist(last_only$rank_at_locus, breaks=0:max(last_only$rank_at_locus))

#  We'd also like to know, out of all genes tested in the region, how close is the gene compared to other nearby ones?
# We obtained this file from BioMart
gene_locs = read.csv("/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/ensembl_genes_with_coordinates.txt", header=TRUE, sep="\t")

gene_locs$Gene.stable.ID = as.character(gene_locs$Gene.stable.ID)
gene_locs$Chromosome.scaffold.name = as.numeric(as.character(gene_locs$Chromosome.scaffold.name))
gene_locs = gene_locs[!is.na(gene_locs$Chromosome.scaffold.name),]

# Get ensembl ID of nearest gene to the SNP in question
priority$nearest = sapply(1:dim(priority)[1], function(x)
       {
		coloc_sub = priority[x,]
		coloc_chrom = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][1])
		coloc_pos = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][2])

       		# Get every nearby gene
       		gene_sub = gene_locs[gene_locs$Chromosome.scaffold.name == coloc_chrom,]
		# Subset down to only the set of genes that we actually tested for colocalization
	        gene_sub = gene_sub[gene_sub$Gene.stable.ID %in% finemap$ensembl,] 

		uniq = gene_sub %>% group_by(Gene.stable.ID) %>% summarize(start = min(c(Gene.start..bp., Gene.end..bp.)), end = max(c(Gene.start..bp., Gene.end..bp.)))

		# For all nearby genes, determine order of proximity to our SNP of interest

		dist = sapply(1:dim(uniq)[1], function(y)
		       {
				this_gene = uniq[y,]
		       		if (this_gene$start > coloc_pos)
				{
					return(this_gene$start - coloc_pos)
				}
				else if (this_gene$end < coloc_pos)
				{
					return(coloc_pos - this_gene$end)
				}
				else
				{
					return(0)
				}
		       }
		)

		uniq = uniq[order(dist),]

		return(uniq$Gene.stable.ID[1])
       })

priority$nearest_hgnc = sapply(priority$nearest, function(x)
	{
		return(as.character(genes[genes$ensembl == x,]$hgnc[1]))
	})
priority$nearest_hgnc[is.na(priority$nearest_hgnc)] = priority$ensembl[is.na(priority$nearest_hgnc)] 
priority$nearest_hgnc[priority$nearest_hgnc == ""] = priority$ensembl[priority$nearest_hgnc == ""] 

# Proximity score = proximity rank of the gene to our SNP, compared
# with all genes tested
priority$proximity = sapply(1:dim(priority)[1], function(x)
       {
		coloc_sub = priority[x,]
		coloc_chrom = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][1])
		coloc_pos = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][2])

       		# Get every nearby gene
       		gene_sub = gene_locs[gene_locs$Chromosome.scaffold.name == coloc_chrom,]
		# Subset down to only the set of genes that we actually tested for colocalization
	        gene_sub = gene_sub[gene_sub$Gene.stable.ID %in% finemap$ensembl,] 

		uniq = gene_sub %>% group_by(Gene.stable.ID) %>% summarize(start = min(c(Gene.start..bp., Gene.end..bp.)), end = max(c(Gene.start..bp., Gene.end..bp.)))

		# For all nearby genes, determine order of proximity to our SNP of interest

		dist = sapply(1:dim(uniq)[1], function(y)
		       {
				this_gene = uniq[y,]
		       		if (this_gene$start > coloc_pos)
				{
					return(this_gene$start - coloc_pos)
				}
				else if (this_gene$end < coloc_pos)
				{
					return(coloc_pos - this_gene$end)
				}
				else
				{
					return(0)
				}
		       }
		)

		uniq = uniq[order(dist),]

		return(which(uniq$Gene.stable.ID == coloc_sub$ensembl))
       })

# Get gene distance from TS
priority$tss_distance = sapply(1:dim(priority)[1], function(x)
       {
		coloc_sub = priority[x,]
		coloc_chrom = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][1])
		coloc_pos = as.numeric(strsplit(coloc_sub$ref_snp, "_")[[1]][2])

		my_gene = gene_locs[gene_locs$Gene.stable.ID == coloc_sub$ensembl,]
		return(min(abs(coloc_pos - my_gene$Transcription.start.site..TSS.)))
       })

# Could ask this question:
# If we consider all the genes that are colocalized with our SNP, what is the proximity
# rank of the BEST one? (This way, we're not penalizing genes that happen to just have a slightly
# lower colocalization score with the nearby gene than with the distant one).

min_proximity = priority %>% group_by(ref_snp) %>% summarize(min_proximity = min(proximity))
hist(min_proximity$min_proximity, breaks=0:max(min_proximity$min_proximity))

#
# Then, for each SNP-gene pair passing the test, show full summary statistics across
# all eQTL files and GWAS
#
gwas_traits = unique(finemap$gwas_trait)

my_matrix = array(0, c(dim(priority)[1], length(gwas_traits)))
#rownames(my_matrix) = paste(priority$ref_snp, " ", priority$nearest_hgnc, " (", priority$hgnc, ", ", priority$proximity, ")", sep="")
rownames(my_matrix) = paste(priority$nearest_hgnc, " --> ", priority$hgnc, sep="")
colnames(my_matrix) = gwas_traits
trunc_matrix = my_matrix	# This matrix will only include the SNPs if they pass the p-value cutoffs
# For each gene that had at least one colocalization
for (i in 1:dim(priority)[1])
{
	# For each GWAS file
	for (k in 1:length(gwas_traits))
	{
		# Fill matrix accordingly, or fill with NA
		result = finemap[finemap$gwas_trait == gwas_traits[k] & 
				 finemap$feature == final_set$feature[i],]
		result_trunc = finemap[finemap$gwas_trait == gwas_traits[k] & 
				 finemap$feature == final_set$feature[i] &
				 finemap$X.log_gwas_pval > 5 & finemap$X.log_eqtl_pval > 5,]

		if (dim(result)[1] == 0)
		{
			my_matrix[i,k] = NA
		}
		else
		{
			my_matrix[i,k] = max(result$clpp)
		}
		if (dim(result_trunc)[1] == 0)
		{
			trunc_matrix[i,k] = NA
		}
		else
		{
			trunc_matrix[i,k] = max(result_trunc$clpp)
		}
	}
}

my_matrix = as.data.frame(my_matrix)

# Clustering by rows and by columns
# Reorder matrix for plotting heatmap
na_status = is.na(my_matrix)
my_matrix[na_status] = 0        # For now, let NAs simply mean no association
hc_hit_snps = hclust(dist(my_matrix))
hc_hit_studies = hclust(dist(t(my_matrix)))
my_matrix[na_status] = NA
my_matrix = my_matrix[hc_hit_snps$order,]
my_matrix = my_matrix[,hc_hit_studies$order]
my_matrix$gene = rownames(my_matrix)

# Plot heatmap
heat = melt(my_matrix)
heat$plot_vals = heat$value
heat$binned_vals = NA
heat$binned_vals[heat$plot_vals < 1.1] = 6
heat$binned_vals[heat$plot_vals < 0.3] = 5
heat$binned_vals[heat$plot_vals < 0.1] = 4
heat$binned_vals[heat$plot_vals < 0.05] = 3
heat$binned_vals[heat$plot_vals < 0.02] = 2
heat$binned_vals[heat$plot_vals < 0.01] = 1
heat$binned_vals = as.factor(heat$binned_vals)

heat$gene = factor(heat$gene, levels=unique(heat$gene))
heat$variable = factor(heat$variable, levels=unique(heat$variable))

heat$strike = ""
heat$strike[is.na(heat$binned_vals)] = "X"
heat$binned_vals[is.na(heat$binned_vals)] = 1

g = ggplot(data = heat, aes(x = variable, y = gene)) +
		geom_tile(aes(fill = binned_vals), color="grey80", size=0.5) +
		scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(6), name="CLPP", labels = c("0-0.01", "0.01-0.02", "0.02-0.05", "0.05-0.1", "0.1-0.3", "0.3-1"),na.value = "white") +
		theme(axis.text.y=element_text(size=10),
		axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(y = "Nearest Gene --> Colocalized Gene") +
		labs(x = "GWAS Trait") + 
		geom_text(label=heat$strike, size=4, color="red", alpha=0.4)

ggsave("/users/mgloud/projects/insulin_resistance/output/insulin_liver_full_clpp.png", width = 8, height = 12, units = "in", dpi = 300, limitsize=FALSE)

trunc_matrix = as.data.frame(trunc_matrix)

# Clustering by rows and by columns
# Reorder matrix for plotting heatmap
# Just cluster them by the full matrix, not the truncated one
trunc_matrix = trunc_matrix[hc_hit_snps$order,]
trunc_matrix = trunc_matrix[,hc_hit_studies$order]
trunc_matrix$gene = rownames(trunc_matrix)

# Plot heatmap
heat = melt(trunc_matrix)
heat$plot_vals = heat$value
heat$binned_vals = NA
heat$binned_vals[heat$plot_vals < 1.1] = 6
heat$binned_vals[heat$plot_vals < 0.3] = 5
heat$binned_vals[heat$plot_vals < 0.1] = 4
heat$binned_vals[heat$plot_vals < 0.05] = 3
heat$binned_vals[heat$plot_vals < 0.02] = 2
heat$binned_vals[heat$plot_vals < 0.01] = 1
heat$binned_vals = as.factor(heat$binned_vals)

heat$gene = factor(heat$gene, levels=unique(heat$gene))
heat$variable = factor(heat$variable, levels=unique(heat$variable))

heat$strike = ""
heat$strike[is.na(heat$binned_vals)] = "X"
heat$binned_vals[is.na(heat$binned_vals)] = 1

g = ggplot(data = heat, aes(x = variable, y = gene)) +
		geom_tile(aes(fill = binned_vals), color="grey80", size=0.5) +
		scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(6), name="CLPP", labels = c("0-0.01", "0.01-0.02", "0.02-0.05", "0.05-0.1", "0.1-0.3", "0.3-1"),na.value = "white") +
		theme(axis.text.y=element_text(size=10),
		axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(y = "Nearest Gene --> Colocalized Gene") +
		labs(x = "GWAS Trait") + 
		geom_text(label=heat$strike, size=4, color="red", alpha=0.4)

ggsave("/users/mgloud/projects/insulin_resistance/output/insulin_liver_with_NAs_clpp.png", width = 8, height = 12, units = "in", dpi = 300, limitsize=FALSE)

stopifnot(FALSE)

# Get total number of genes tested at this locus
genes_tested = finemap %>% group_by(rsid) %>% summarize(genes_tested=length(unique(as.character(ensembl))))
priority = merge(priority, genes_tested)

# Get total number of colocalizations found at this locus
last_only$num_colocs = last_only$rank_at_locus
priority = merge(priority, last_only[c("rsid", "num_colocs")], by=c("rsid"))

# Get ALL hits at this locus (regardless of whether others were found in same tissue
extra_info = priority[c("rsid", "ensembl", "rank_at_locus", "proximity", "tss_distance", "num_colocs", "genes_tested")]
priority = merge(pre_final_set, extra_info, by=c("rsid", "ensembl"))

num_gwas_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_gwas_colocalized = length(unique(base_gwas_file)))
priority = merge(priority, num_gwas_colocalized, by=c("rsid", "ensembl"))

num_eqtl_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_tissues_colocalized = length(unique(eqtl_file)))
priority = merge(priority, num_eqtl_colocalized, by=c("rsid", "ensembl"))

num_gwas_by_eqtl_colocalized = priority %>% group_by(rsid, ensembl) %>% summarize(num_gwas_by_eqtl_colocalized = length(unique(paste(base_gwas_file, eqtl_file))))
priority = merge(priority, num_gwas_by_eqtl_colocalized, by=c("rsid", "ensembl"))

priority = priority[rev(order(priority$clpp_mod)),]
write.table(priority, file="/users/mgloud/projects/brain_gwas/scripts/auxiliary/eriks_grant/results/prioritized_gene_rankings_no_clpp_mod_cutoff.tsv", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
