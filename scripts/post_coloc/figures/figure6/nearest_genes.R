require(dplyr)
require(readr)

# Code to get the nearest genes for the column in Figure 6

results = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", sep="\t", header=TRUE)

# I got rid of this for now because I changed it to be for hg38 instead. Really though, I should adapt this
# code to work for both
ensembl_locs = read_delim("data/gencode/gencode.v31.annotation.gtf.gz", delim="\t", col_names=FALSE, skip=5)
ensembl_locs = ensembl_locs[,c(1,3,4,5,7,9)]
colnames(ensembl_locs) = c("chrom", "type", "gene_start", "gene_end", "strand", "prefeature")

tmp = ensembl_locs[ensembl_locs$strand == "-",]$gene_start
ensembl_locs[ensembl_locs$strand == "-",]$gene_start = ensembl_locs[ensembl_locs$strand == "-",]$gene_end
ensembl_locs[ensembl_locs$strand == "-",]$gene_end = tmp

ensembl_locs$feature = sapply(ensembl_locs$prefeature, function(x) {strsplit(x, "gene_id\ \"")[[1]][2]})
ensembl_locs$feature = sapply(ensembl_locs$feature, function(x) {strsplit(x, "\\.")[[1]][1]})
ensembl_locs$hgnc = sapply(ensembl_locs$prefeature, function(x) {strsplit(x, "gene_name\ \"")[[1]][2]})
ensembl_locs$hgnc = sapply(ensembl_locs$hgnc, function(x) {strsplit(x, "\"")[[1]][1]})
ensembl_locs = ensembl_locs[ensembl_locs$type == "gene",]

results = results %>% filter(clpp_mod > 0.35)

# Get distance from coloc SNPs to TSS's of coloc'ed genes
all_snps = results %>% select(chr, pos) %>% filter(!duplicated(paste(chr,pos)))

deets = sapply(1:dim(all_snps)[1], function(x)
{
		answers = list()

	       row = all_snps[x,]
	       top_left_genes = ensembl_locs %>% filter(chrom == paste0("chr", row$chr)) %>% filter(pmax(row$pos - gene_start, row$pos - gene_end) > 0) %>% mutate(left_dist=pmax(0,pmin(row$pos - gene_start, row$pos - gene_end))) %>% arrange(left_dist)
	       top_right_genes = ensembl_locs %>% filter(chrom == paste0("chr", row$chr)) %>% filter(pmax(gene_start - row$pos, gene_end - row$pos) > 0) %>% mutate(right_dist=pmax(0,pmin(gene_start - row$pos, gene_end - row$pos))) %>% arrange(right_dist)

	       answers$left_dist = top_left_genes$left_dist[1]
	       answers$left_gene = top_left_genes$hgnc[1]
	       answers$left_feature = top_left_genes$feature[1]
	       answers$right_dist = top_right_genes$right_dist[1]
	       answers$right_gene = top_right_genes$hgnc[1]
	       answers$right_feature = top_right_genes$feature[1]

	       return(answers)
       }
)

deets = as.data.frame(t(deets))
deets$left_dist = unlist(deets$left_dist)
deets$left_gene = unlist(deets$left_gene)
deets$left_feature = unlist(deets$left_feature)
deets$right_dist = unlist(deets$right_dist)
deets$right_gene = unlist(deets$right_gene)
deets$right_feature = unlist(deets$right_feature)

all_snps = cbind(all_snps, deets)

write.table(all_snps, file="output/post_coloc/plots/figure6/nearest_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

