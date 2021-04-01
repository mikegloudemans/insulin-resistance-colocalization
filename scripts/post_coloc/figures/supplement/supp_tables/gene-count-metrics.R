# Counts of genes / loci tested and colocalized for a few of the supplementary tables

clpp_mod_cutoff = 0.35

d = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", header=TRUE, sep="\t")

unique(d$gwas_short)
sapply(unique(d$gwas_short), function(x) {length(unique(d[d$gwas_short == x,]$hgnc))})
sapply(unique(d$gwas_short), function(x) {length(unique(d[d$gwas_short == x,]$locus))})
sapply(unique(d$gwas_short), function(x) {length(unique(d[d$gwas_short == x & d$clpp_mod >= clpp_mod_cutoff,]$hgnc))})
sapply(unique(d$gwas_short), function(x) {length(unique(d[d$gwas_short == x & d$clpp_mod >= clpp_mod_cutoff,]$locus))})

unique(d$eqtl_short)
sapply(unique(d$eqtl_short), function(x) {length(unique(d[d$eqtl_short == x,]$hgnc))})
sapply(unique(d$eqtl_short), function(x) {length(unique(d[d$eqtl_short == x,]$locus))})
sapply(unique(d$eqtl_short), function(x) {length(unique(d[d$eqtl_short == x & d$clpp_mod >= clpp_mod_cutoff,]$hgnc))})
sapply(unique(d$eqtl_short), function(x) {length(unique(d[d$eqtl_short == x & d$clpp_mod >= clpp_mod_cutoff,]$locus))})

tissues = c("Subcutaneous-Adipose", "Visceral-Adipose", "Liver", "Pancreas", "Skeletal-Muscle")
sapply(tissues, function(x) {length(unique(d[grepl(x, d$eqtl_short),]$hgnc))})
sapply(tissues, function(x) {length(unique(d[grepl(x, d$eqtl_short),]$locus))})
sapply(tissues, function(x) {length(unique(d[grepl(x, d$eqtl_short) & d$clpp_mod >= clpp_mod_cutoff,]$hgnc))})
sapply(tissues, function(x) {length(unique(d[grepl(x, d$eqtl_short) & d$clpp_mod >= clpp_mod_cutoff,]$locus))})

