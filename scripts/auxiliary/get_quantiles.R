## Get quantiles of CLPP and CLPP-mod scores, so that
## we have a better idea of what a given cutoff means

require(qvalue)

coloc_results_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/full_coloc_results_qced.txt" 

# TODO: Fix path
results = read.table(coloc_results_file, sep="\t", header=TRUE)
results$gwas_p = 10^(-results$neg_log_gwas_pval)
results$eqtl_p = 10^(-results$neg_log_eqtl_pval)

sqtl = results[grepl("sQTL", results$eqtl_file),]
eqtl = results[grepl("eQTL", results$eqtl_file),]


# Quantiles for various cutoffs
quantile(results$clpp_mod, seq(0,1,by=0.05))
quantile(sqtl$clpp_mod, seq(0,1,by=0.05))
quantile(eqtl$clpp_mod, seq(0,1,by=0.05))

# Histogram of values
hist(1-results$clpp)
hist(1-results$clpp_mod)


