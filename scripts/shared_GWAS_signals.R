# Author: Mike Gloudemans
# Date: 3/9/2018

require(readr)
require(ggplot2)
require(reshape2)
require(data.table)
require(ggdendro)
require(cowplot)

# Noticeably, "HbA1c_MAGIC_Mixed_AllSNPs.txt" is not included in this
# list of summary stats. That's because its significance values
# were listed as Bayes factors, and we didn't figure out yet how to
# turn these into p-values most effectively.

# I've intentionally excluded MI traits from these summary statistics
# because the stats for these aren't true summary statistics
sum_stats_files = c("2hrGlu_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
		    "BMI_GIANT_Europeans_AllSNPs.prepared.txt.gz",
		    "CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt.gz",
		    "BMI_GIANT_Mixed_AllSNPs.prepared.txt.gz",
		    "FastGlu_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
		    "FastInsu_adjBMI_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
		    "FastInsu_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
		    "HDL_GCLC_Mixed.prepared.txt.gz",
		    "IS_METASTROKE_Mixed.prepared.txt.gz",
		    "MAGIC_ISI_Model_1_AgeSexOnly.prepared.txt.gz",
		    "MAGIC_ISI_Model_2_AgeSexBMI.prepared.txt.gz",
		    "T2D_DIAGRAM_European.prepared.txt.gz",
		    "TG_GCLC_Mixed.prepared.txt.gz",
		    "WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt.gz",
		    "WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt.gz"
		    )

# I'm not sure why we have "MAGIC_ISI_Model_1" and "MAGIC_ISI_Model_2" GWAS,
# but we don't have top hits ISI1. Probably not a huge concern though.
# Also, I commented out the files for which we have top hits but no summary statistics,
# because they mostly just clutter up the plot with empty space.
hits_files = c("/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/CHD_CARDIoGRAMplusC4D_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/2hrGlu_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/IS_METASTROKE_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/FastInsu_adjBMI_MAGIC_Europeans_Hits_hg19.txt",
	       #"/users/mgloud/projects/insulin_resistance/data/MI_adjBMI_GENESISGuardian_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/FastInsu_MAGIC_Europeans_Hits_hg19.txt",
	       #"/users/mgloud/projects/insulin_resistance/data/HbA1c_MAGIC_Mixed_Hits_hg19.txt",
	       #"/users/mgloud/projects/insulin_resistance/data/MI_GENESISGuardian_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/HDL_GCLC_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/TG_GCLC_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/FastGlu_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/ISI2_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/T2D_DIAGRAM_European_Hits_hg19.txt"
	       )

window = 100000


sum_stats = list()

for (ss in sum_stats_files)
{
	print(ss)
	sum_stats[[ss]] = read_delim(gzfile(paste("/users/mgloud/projects/gwas/data/prepared/", ss, sep="")), delim="\t")
	sum_stats[[ss]]$chr = gsub("chr", "", sum_stats[[ss]]$chr)
}


# TODO: Create matrix for storing GWAS results across traits
# Actually create two
hit_pvals = array(rep(0, 10000 * length(names(sum_stats))), dim=c(10000, length(names(sum_stats))))
window_pvals = array(rep(0, 10000 * length(names(sum_stats))), dim=c(10000, length(names(sum_stats))))

rownames(hit_pvals) = 1:10000
rownames(window_pvals) = 1:10000

colnames(hit_pvals) = sapply(strsplit(names(sum_stats), "\\."), function(x){x[[1]][1]})
colnames(window_pvals) = sapply(strsplit(names(sum_stats), "\\."), function(x){x[[1]][1]})

hit_num = 1
for (GWAS in hits_files)
{
	print(GWAS)
	hits = read_delim(GWAS, delim="\t")
	hits$chr = gsub("chr", "", hits$chr)

	if (dim(hits)[1] == 0)
	{
		next
	}

	for (i in 1:dim(hits)[1])
	{
		print(i)
		# Loop through all summary stats files

		for (j in 1:length(names(sum_stats)))
		{
			ss = names(sum_stats)[j]
			other_GWAS = sum_stats[[ss]]

			chr_sub = other_GWAS[other_GWAS$chr == hits$chr[i],]

			SNP = chr_sub[chr_sub$snp_pos == hits$snp_pos[i],]
			if (length(SNP$pvalue) == 0)
			{
				hit_pvalue = NA
			}
			else
			{
				hit_pvalue = min(SNP$pvalue)
			}

			SNPs = chr_sub[chr_sub$snp_pos > (hits$snp_pos[i] - window) & 
					 chr_sub$snp_pos < (hits$snp_pos[i] + window),]
			if (length(SNPs$pvalue) == 0)
			{
				window_pvalue = NA
			}
			else
			{
				window_pvalue = min(SNPs$pvalue)
			}

			hit_pvals[hit_num, j] = hit_pvalue
			window_pvals[hit_num, j] = window_pvalue

		}

		trait = strsplit(GWAS, "/")[[1]]
		trait = trait[length(trait)]
		trait = strsplit(trait, "\\.")
		trait = trait[[1]][1]
		rownames(hit_pvals)[hit_num] = paste(hits$chr[i], hits$snp_pos[i], trait, sep="_")
		rownames(window_pvals)[hit_num] = paste(hits$chr[i], hits$snp_pos[i], trait, sep="_")
		
		hit_num = hit_num + 1
	}
}

hit_pvals = hit_pvals[1:(hit_num-1),]
window_pvals = window_pvals[1:(hit_num-1),]

stuff = list(hit_pvals, window_pvals)
save(stuff, file="myfile.RData")

load("myfile.RData")
hit_pvals = stuff[[1]]
window_pvals = stuff[[2]]

# A possible way to remove redundancy:
# Identify SNPs in heatmap only by position; not by study.
# If two SNPs are within 10000 bp of one another, remove one of them.

loci = as.data.frame(list(chrom = sapply(rownames(hit_pvals), function(x) {strsplit(x, "_")[[1]][1]}), 
			  pos = sapply(rownames(hit_pvals), function(x) {as.numeric(strsplit(x, "_")[[1]][2])})))
# Return a boolean value stating whether this SNP is a duplicate of one above it in the list
dup_snps = sapply(1:dim(loci)[1], function(x)
	{
		my_chrom = loci$chrom[x]
		my_pos = loci$pos[x]

		return(sum((loci$chrom[1:x] == my_chrom) & (abs(loci$pos[1:x] - my_pos) < 100000)) > 1)
	})

# If SNP isn't on a conventional chromosome, just get rid of it
dup_snps[is.na(dup_snps)] = TRUE

hit_pvals = hit_pvals[!dup_snps,]
window_pvals = window_pvals[!dup_snps,]
loci = loci[!dup_snps,]

# I don't think any of the summary stats I was given have the X chromosome,
# so just remove those SNPs for now
is_x = loci$chr == "X"
hit_pvals = hit_pvals[!is_x,]
window_pvals = window_pvals[!is_x,]
loci = loci[!is_x,]

# Rename SNPs so that they're not specific to any one single study
rownames(hit_pvals) = paste(loci$chrom, loci$pos, sep="_")
rownames(window_pvals) = paste(loci$chrom, loci$pos, sep="_")

# Note: This is just for visualization...comparing pvalues from one
# study to another isn't a rigorous method of comparison, and I'll
# have to think carefully about how now to make this plot misleading.
hit_pvals = -log10(hit_pvals)
hit_pvals = pmin(hit_pvals, 20)	# At a certain point, a hit is a hit 
window_pvals = -log10(window_pvals)
window_pvals = pmin(window_pvals, 20)	# At a certain point, a hit is a hit 

# Now, get rid of all SNPs for which p-vals aren't at least 
# borderline significant in even a single study.
window_pvals = window_pvals[apply(window_pvals, 1, max, na.rm=TRUE) > 5,]

# Clustering by rows and by columns
# Reorder matrix for plotting heatmap
na_status = is.na(hit_pvals)
hit_pvals[na_status] = 0	# For now, let NAs simply mean no association
hc_hit_snps = hclust(dist(hit_pvals))
hc_hit_studies = hclust(dist(t(hit_pvals)))
hit_pvals[na_status] = NA
hit_pvals = hit_pvals[hc_hit_snps$order,]
hit_pvals = hit_pvals[,hc_hit_studies$order]

na_status = is.na(window_pvals)
window_pvals[na_status] = 0	# For now, let NAs simply mean no association
hc_window_snps = hclust(dist(window_pvals))
hc_window_studies = hclust(dist(t(window_pvals)))
window_pvals[na_status] = NA
window_pvals = window_pvals[hc_window_snps$order,]
window_pvals = window_pvals[,hc_window_studies$order]

# Plot dendrograms. Later, we'll assemble the plot by code; for now, we'll just have
# to stitch them together manually
g_hsn_dendro = ggdendrogram(hc_hit_snps) + coord_flip() + scale_y_reverse() + theme_dendro()
g_hsn_dendro
#ggsave("gwas_tophit_snp_dendro.png", width = 8, height = 100, units = "in", dpi = 300, limitsize=FALSE)

g_hst_dendro = ggdendrogram(hc_hit_studies) + scale_y_reverse() + theme_dendro()
g_hst_dendro
#ggsave("gwas_tophit_study_dendro.png", width = 8, height = 8, units = "in", dpi = 300, limitsize=FALSE)

g_wsn_dendro = ggdendrogram(hc_window_snps) + coord_flip() + scale_y_reverse() + theme_dendro()
g_wsn_dendro
#ggsave("gwas_window_snp_dendro.png", width = 8, height = 100, units = "in", dpi = 300, limitsize=FALSE)

g_wst_dendro = ggdendrogram(hc_window_studies) + scale_y_reverse() + theme_dendro()
g_wst_dendro
#ggsave("gwas_window_study_dendro.png", width = 8, height = 8, units = "in", dpi = 300, limitsize=FALSE)



# Plot heatmap
heat = melt(hit_pvals, id="gene")
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_text(size=7),
		      axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(x = "Replication GWAS") +
		labs(y = "Discovery GWAS Hit")

ggsave("gwas_tophit_replication.png", width = 8, height = 50, units = "in", dpi = 300, limitsize=FALSE)

# Y-axis labels removed and heatmap highly compressed
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_blank(),
		      axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(x = "Replication GWAS") +
		labs(y = "Discovery GWAS Hit")

ggsave("gwas_tophit_replication_compressed.png", width = 4, height = 10, units = "in", dpi = 300, limitsize=FALSE)

# Same thing for window-adjusted
heat = melt(window_pvals, id="gene")
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_text(size=7),
		      axis.text.x = element_text(angle = 90, hjust = 1)) +
		labs(x = "Replication GWAS") +
		labs(y = "Discovery GWAS Hit")

ggsave("gwas_window_replication.png", width = 8, height = 50, units = "in", dpi = 300, limitsize=FALSE)

# Y-axis labels removed and heatmap highly compressed
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_blank(),
		      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
		      axis.title.x=element_blank(),
		      axis.title.y=element_blank())

ggsave("gwas_window_replication_compressed.png", width = 4, height = 10, units = "in", dpi = 300, limitsize=FALSE)

# compressed windows, with dendrograms

combined = ggdraw() + 
		draw_plot(g, 0.4, 0.2, 0.6, 0.8) +
		draw_plot(g_wsn_dendro, 0, 0.377, 0.42, 0.65) +
		draw_plot(g_wst_dendro, 0.395, 0, 0.51, 0.215)

ggsave("gwas_window_replication_compressed_combined.png", width = 8, height = 16, units = "in", dpi = 300, limitsize=FALSE)



