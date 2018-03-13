# Author: Mike Gloudemans
# Date: 2/2/2018

require(readr)
require(ggplot2)
require(reshape2)

main = function()
{
	hits_files = c("/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Mixed_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Europeans_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Mixed_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/CHD_CARDIoGRAMplusC4D_Mixed_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/IS_METASTROKE_Mixed_Hits_hg19.txt",
		       "/users/mgloud/projects/insulin_resistance/data/T2D_DIAGRAM_European_Hits_hg19.txt"
		       )
	
	finemap_files = c("/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/BMI_GIANT_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/BMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/WHRadjBMI_GIANT_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/WHRadjBMI_GIANT_Mixed_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/CHD_CARDIoGRAMplusC4D_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/IS_METASTROKE_Mixed_prepared_txt_gz_finemap_clpp_status.txt",
			  "/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/T2D_DIAGRAM_European_prepared_txt_gz_finemap_clpp_status.txt"

			  #"/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/FastGlu_MAGIC_Europeans_AllSNPs_prepared_txt_gz_finemap_clpp_status.txt",
			  #"/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/FastInsu_adjBMI_MAGIC_Europeans_Hits_hg19.txt",
			  #"/users/mgloud/projects/brain_gwas/output/completed/insulin/2018-01-28_08-55-30_insulin_gtex_for_bruna/FastInsu_MAGIC_Europeans_Hits_hg19.txt"
		       )

	out_files = c("/users/mgloud/projects/insulin_resistance/output/plots/BMI_GIANT_Europeans_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/BMI_GIANT_Mixed_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/WHRadjBMI_GIANT_Europeans_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/WHRadjBMI_GIANT_Mixed_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/CHD_CARDIoGRAMplusC4D_Mixed_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/IS_METASTROKE_Mixed_Hits_clpp_heatmap.png",
		      "/users/mgloud/projects/insulin_resistance/output/plots/T2D_DIAGRAM_European_Hits_clpp_heatmap.png"
		      )

	input_files = data.frame(list(hits=hits_files, finemap=finemap_files, outs = out_files), stringsAsFactors=FALSE)

	tissues = c("Adipose_Visceral_Omentum", "Adipose_Subcutaneous", "Artery_Aorta", "Liver", "Artery_Coronary", "Muscle_Skeletal")

	for (i in 1:dim(input_files)[1])
	{
		print(i)
		plot_clpp_heatmap(input_files$hits[i], input_files$finemap[i], input_files$outs[i], tissues)
	}

}

plot_clpp_heatmap = function(hits_file, finemap_file, out_file, tissues)
{
	print(hits_file)
	hits = read_delim(hits_file, delim="\t")
	hits$chr = as.numeric(gsub("chr", "", hits$chr))
	hits = hits[!is.na(hits$chr),]

	print(finemap_file)
	finemap = read_delim(finemap_file, delim="\t", col_names=c("snp", "eqtl_file", "gwas_file", "gene", "conditional_level", "n_snps", "clpp", "neg_logp_gwas", "neg_logp_eqtl"))
	finemap$chrom = sapply(finemap$snp, function(x){as.numeric(strsplit(x, "_")[[1]][1])})
	finemap$pos = sapply(finemap$snp, function(x){as.numeric(strsplit(x, "_")[[1]][2])})

	# Since we don't know the target size, just make it larger than needed
	all_coloc_scores = array(0, c(1000, length(tissues)))
	colnames(all_coloc_scores) = tissues
	rownames(all_coloc_scores) = 1:1000
	coloc_num = 0
	for (i in 1:dim(hits)[1])
	{
		print(i)

		hit = hits[i,]

		hit_coloc = finemap[finemap$chrom == hit$chr & abs(finemap$pos - hit$snp_pos) <= 500000,]
		hit_coloc = hit_coloc[hit_coloc$neg_logp_gwas > 5 & hit_coloc$neg_logp_eqtl > 5,]

		if (dim(hit_coloc)[1] == 0)
		{ 	
			coloc_num = coloc_num + 1
			rownames(all_coloc_scores)[coloc_num] = paste(hit$chr, hit$snp_pos, sep="_")
			all_coloc_scores[coloc_num,] = rep(NA, length(tissues))
			next
		}
		print(max(hit_coloc$clpp))	

		# Figure out which genes have putative colocalization in any tissue
		best_genes = unique(hit_coloc$gene[hit_coloc$clpp > 0.01])
		if (length(best_genes) == 0)
		{
			best_genes = hit_coloc$gene[which.max(hit_coloc$clpp)]
		}
		print(best_genes)

		for (bg in best_genes)
		{
			# Skip this SNP/gene combo if we've already analyzed it for a different hit...
			if (paste(hit$chr, hit$snp_pos, bg, sep="_") %in% rownames(all_coloc_scores))
			{
				next
			}
			coloc_num = coloc_num + 1
			rownames(all_coloc_scores)[coloc_num] = paste(hit$chr, hit$snp_pos, bg, sep="_")	

			# Get an unfiltered version of the hits so that we we can see scores for 
			# genes that are just barely sub-threshold
			hit_unfiltered_coloc = finemap[finemap$chrom == hit$chr & abs(finemap$pos - hit$snp_pos) <= 500000,]
			gene_coloc = hit_unfiltered_coloc[hit_unfiltered_coloc$gene == bg,]

			# See to what extent that gene colocalizes in each tissue, and store these
			# results in a matrix for heatmap plotting.
			
			for (j in 1:length(tissues))
			{
				tissue = paste(tissues[j], "_allpairs_txt_gz", sep="")
				tissue_coloc = gene_coloc[gene_coloc$eqtl_file == tissue,]
				if (dim(tissue_coloc)[1] == 0)
				{
					all_coloc_scores[coloc_num,j] = NA
				}
				else
				{
					# Sanity check: we should only have one score left
					# at this point
					stopifnot(dim(tissue_coloc)[1] == 1)
					all_coloc_scores[coloc_num,j] = tissue_coloc$clpp
				}
			}
		}
	}

	# Discard extra allocated rows
	all_coloc_scores = all_coloc_scores[1:(coloc_num-1),]

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
			scale_fill_manual(drop=FALSE, values=colorRampPalette(c("white","red"))(5), name="CLPP", labels = c("0-0.005", "0.005-0.01", "0.01-0.02", "0.02-0.1", "0.1-1"),na.value = "grey80") +
			theme(axis.text.y=element_text(size=7),
			      axis.text.x = element_text(angle = 45, hjust = 1)) +
			labs(x = "Tissue") +
			labs(y = "SNP")

	g

	ggsave(out_file, width = 4, height = (1/8) * dim(fix)[1], units = "in", dpi = 300)
}


main()

# TODO: See if there are any other genes colocalizing at the same loci,
# but overshadowed by top variant

# TODO: See how many of the loci in the full coloc results file that show
# high probabilities of colocalization, aren't listed in this list due to
# not being near top genes

# TODO: Aggregate colocalization across all hits in a tissue to see if certain tissues
# are enriched for colocalization
