require(dplyr)
require(readr)

results = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", sep="\t", header=TRUE)

# I got rid of this for now because I changed it to be for hg38 instead. Really though, I should adapt this
# code to work for both
ensembl_locs = read_delim("data/gencode/gencode.v31.annotation.gtf.gz", delim="\t", col_names=FALSE, skip=5)
ensembl_locs = ensembl_locs[,c(1,3,4,5,7,9)]
colnames(ensembl_locs) = c("chrom", "type", "gene_start", "gene_end", "strand", "feature")

tmp = ensembl_locs[ensembl_locs$strand == "-",]$gene_start
ensembl_locs[ensembl_locs$strand == "-",]$gene_start = ensembl_locs[ensembl_locs$strand == "-",]$gene_end
ensembl_locs[ensembl_locs$strand == "-",]$gene_end = tmp

ensembl_locs$feature = sapply(ensembl_locs$feature, function(x) {strsplit(x, "gene_id\ \"")[[1]][2]})
ensembl_locs$feature = sapply(ensembl_locs$feature, function(x) {strsplit(x, "\\.")[[1]][1]})
ensembl_locs = ensembl_locs[ensembl_locs$type == "gene",]

results$gene_start = ensembl_locs[match(results$ensembl, ensembl_locs$feature),]$gene_start
results$gene_end = ensembl_locs[match(results$ensembl, ensembl_locs$feature),]$gene_end
results$strand = ensembl_locs[match(results$ensembl, ensembl_locs$feature),]$strand

results = results %>% filter(clpp_mod > 0.35)

# Get distance from coloc SNPs to TSS's of coloc'ed genes
all_tss_dists = results %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))
eqtl_tss_dists = results %>% filter(clpp_mod > 0.35) %>% filter(grepl("eQTL", eqtl_file)) %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))
sqtl_tss_dists = results %>% filter(clpp_mod > 0.35) %>% filter(grepl("sQTL", eqtl_file)) %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))

print(head(all_tss_dists))

all_single_coloc_tss_dists = results %>% filter(clpp_mod > 0.35) %>% filter((step1 == "loci_2") | (step1 == "loci_0.1")) %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))
eqtl_single_coloc_tss_dists = results %>% filter(clpp_mod > 0.35) %>% filter((step1 == "loci_2") | (step1 == "loci_0.1")) %>% filter(grepl("eQTL", eqtl_file)) %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))
sqtl_single_coloc_tss_dists = results %>% filter(clpp_mod > 0.35) %>% filter((step1 == "loci_2") | (step1 == "loci_0.1")) %>% filter(grepl("sQTL", eqtl_file)) %>% group_by(chr, pos, ensembl) %>% summarize(tss_dist = (pos[1]-gene_start[1]) * ifelse(strand[1]=="-",-1,1))
breaks = seq(-1e+06, 1e+06, by = 20000)
hist(all_tss_dists$tss_dist, breaks=breaks)
print(median(all_tss_dists$tss_dist, na.rm=TRUE))
abline(v=median(all_tss_dists$tss_dist, na.rm=TRUE), lwd=4, col="red", lty=3)
hist(eqtl_tss_dists$tss_dist, breaks=breaks)
abline(v=median(eqtl_tss_dists$tss_dist, na.rm=TRUE), lwd=4, col="red", lty=3)
hist(sqtl_tss_dists$tss_dist, breaks=breaks)
abline(v=median(sqtl_tss_dists$tss_dist, na.rm=TRUE), lwd=4, col="red", lty=3)
hist(all_single_coloc_tss_dists$tss_dist, breaks=breaks)
abline(v=median(all_single_coloc_tss_dists$tss_dist), lwd=4, col="red", lty=3)
hist(eqtl_single_coloc_tss_dists$tss_dist, breaks=breaks)
abline(v=median(eqtl_single_coloc_tss_dists$tss_dist), lwd=4, col="red", lty=3)
hist(sqtl_single_coloc_tss_dists$tss_dist, breaks=breaks)
abline(v=median(sqtl_single_coloc_tss_dists$tss_dist), lwd=4, col="red", lty=3)

all_tss_dists$coloc_proximity_rank = sapply(1:dim(all_tss_dists)[1], function(x)
       {
	       row = all_tss_dists[x,]
	       top_genes = ensembl_locs %>% filter(chrom == paste0("chr", row$chr)) %>% arrange(pmin(abs(gene_start - row$pos), abs(gene_end - row$pos)))
	       answer = which(top_genes$feature == row$ensembl)
	       if (length(answer) == 0)
	       {
		       answer = NA
	       }
	       return(answer)
       }
)
all_tss_dists = all_tss_dists %>% filter(!is.na(coloc_proximity_rank))
all_tss_dists[all_tss_dists$coloc_proximity_rank > 50,]$coloc_proximity_rank = 50

all_single_coloc_tss_dists$coloc_proximity_rank = sapply(1:dim(all_single_coloc_tss_dists)[1], function(x)
       {
	       row = all_single_coloc_tss_dists[x,]
	       top_genes = ensembl_locs %>% filter(chrom == paste0("chr", row$chr)) %>% arrange(pmin(abs(gene_start - row$pos), abs(gene_end - row$pos)))
	       answer = which(top_genes$feature == row$ensembl)
	       if (length(answer) == 0)
	       {
		       answer = NA
	       }
	       return(answer)
       }
)
all_single_coloc_tss_dists = all_single_coloc_tss_dists %>% filter(!is.na(coloc_proximity_rank))
all_single_coloc_tss_dists[all_single_coloc_tss_dists$coloc_proximity_rank > 50,]$coloc_proximity_rank = 50

hist(all_tss_dists$coloc_proximity_rank, breaks = seq(-0.5,50.5))
hist(all_single_coloc_tss_dists$coloc_proximity_rank, breaks = seq(-0.5,50.5))

