# Author: Mike Gloudemans
# Date: 3/9/2018

require(readr)
require(ggplot2)
require(reshape2)
require(data.table)

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
hits_files = c("/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/WHRadjBMI_GIANT_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/CHD_CARDIoGRAMplusC4D_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/2hrGlu_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/IS_METASTROKE_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/FastInsu_adjBMI_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/MI_adjBMI_GENESISGuardian_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/FastInsu_MAGIC_Europeans_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/HbA1c_MAGIC_Mixed_Hits_hg19.txt",
	       "/users/mgloud/projects/insulin_resistance/data/MI_GENESISGuardian_Mixed_Hits_hg19.txt",
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

# Note: This is just for visualization...comparing pvalues from one
# study to another isn't a rigorous method of comparison, and I'll
# have to think carefully about how now to make this plot misleading.
hit_pvals = -log10(hit_pvals)
hit_pvals = pmin(hit_pvals, 20)	# At a certain point, a hit is a hit 
window_pvals = -log10(window_pvals)
window_pvals = pmin(window_pvals, 20)	# At a certain point, a hit is a hit 

# Plot heatmap
heat = melt(hit_pvals, id="gene")
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_text(size=7),
		      axis.text.x = element_text(angle = 45, hjust = 1)) +
		labs(x = "Replication GWAS") +
		labs(y = "Discovery GWAS Hit")

ggsave("gwas_tophit_replication.png", width = 8, height = 100, units = "in", dpi = 300, limitsize=FALSE)

# Same thing for window-adjusted
heat = melt(window_pvals, id="gene")
g = ggplot(data = heat, aes(x = Var2, y = Var1)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradient2(low = "white", high = 'orangered4', midpoint = median(heat$value, na.rm=TRUE)) +
		theme(axis.text.y=element_text(size=7),
		      axis.text.x = element_text(angle = 45, hjust = 1)) +
		labs(x = "Replication GWAS") +
		labs(y = "Discovery GWAS Hit")

ggsave("gwas_window_replication.png", width = 8, height = 100, units = "in", dpi = 300, limitsize=FALSE)
