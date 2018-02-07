# Author: Mike Gloudemans
# Date: 2/2/2018

require(readr)
require(ggplot2)
require(reshape2)


tissues = c("Adipose_Visceral_Omentum", "Adipose_Subcutaneous", "Artery_Aorta", "Liver", "Artery_Coronary", "Muscle_Skeletal")

hits = read_delim("/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt", delim="\t")
hits$chr = as.numeric(gsub("chr", "", hits$chr))

finemap = read_delim("/users/mgloud/projects/brain_gwas/output/2018-01-28_08-55-30_insulin_gtex_for_bruna/BMI_GIANT_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt", delim="\t", col_names=c("snp", "eqtl_file", "gwas_file", "gene", "conditional_level", "n_snps", "clpp", "neg_logp_gwas", "neg_logp_eqtl"))
finemap$chrom = sapply(finemap$snp, function(x){as.numeric(strsplit(x, "_")[[1]][1])})
finemap$pos = sapply(finemap$snp, function(x){as.numeric(strsplit(x, "_")[[1]][2])})


all_coloc_scores = array(0, c(dim(hits)[1], length(tissues)))
colnames(all_coloc_scores) = tissues
rownames(all_coloc_scores) = paste(hits$chr, "_", hits$snp_pos, sep="")
for (i in 1:dim(hits)[1])
{
	hit = hits[i,]
	hit_coloc = finemap[finemap$chrom == hit$chr & abs(finemap$pos - hit$snp_pos) <= 500000,]
	hit_coloc = hit_coloc[hit_coloc$neg_logp_gwas > 5 & hit_coloc$neg_logp_eqtl > 5,]

	print(i)

	if (dim(hit_coloc)[1] == 0)
	{
		next
	}
	print(max(hit_coloc$clpp))	

	# Figure out which gene has the strongest colocalization in any tissue
	best_gene = hit_coloc$gene[which.max(hit_coloc$clpp)]
	print(best_gene)

	rownames(all_coloc_scores)[i] = paste(rownames(all_coloc_scores)[i], "_", best_gene, sep="")
	

	gene_coloc = hit_coloc[hit_coloc$gene == best_gene,]

	# See to what extent that gene colocalizes in each tissue, and store these
	# results in a matrix for heatmap plotting.
	
	for (j in 1:length(tissues))
	{
		tissue = paste(tissues[j], "_allpairs_txt_gz", sep="")
		tissue_coloc = gene_coloc[gene_coloc$eqtl_file == tissue,]
		if (dim(tissue_coloc)[1] == 0)
		{
			all_coloc_scores[i,j] = 0
		}
		else
		{
			# Sanity check: we should only have one score left
			# at this point
			stopifnot(dim(tissue_coloc)[1] == 1)
			all_coloc_scores[i,j] = tissue_coloc$clpp
		}
	}
}

# Bin scores to make heatmap look nicer
binned_coloc_scores = cut(all_coloc_scores, breaks=c(-0.001, 0.005, 0.01, 0.02, 0.1, 1), labels=1:5)
dim(binned_coloc_scores) = dim(all_coloc_scores)

# Reshape it back into proper dimensions with proper labels and everything
# It's ridiculous that we have to do it this way, but we
# must do what we must do
fix = data.frame(matrix(0, nrow=dim(binned_coloc_scores)[1]))
for (c in colnames(all_coloc_scores))
{
	fix[c] = binned_coloc_scores[,which(colnames(all_coloc_scores) == c)]
}
fix = fix[,-1]
rownames(fix) = rownames(all_coloc_scores)
fix$gene = rownames(fix)

# Plot heatmap
heat = melt(fix, id="gene")
g = ggplot(data = heat, aes(x = variable, y = gene)) +
		geom_tile(aes(fill = value)) +
		scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(5), name="CLPP", labels = c("0-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.1", "0.1-1")) +
		theme(axis.text.y=element_text(size=7),
		      axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(x = "Tissue") +
		labs(y = "SNP")

g

# TODO: See if there are any other genes colocalizing at the same loci,
# but overshadowed by top variant

# TODO: See how many of the loci in the full coloc results file that show
# high probabilities of colocalization, aren't listed in this list due to
# not being near top genes

# TODO: Aggregate colocalization across all hits in a tissue to see if certain tissues
# are enriched for colocalization
