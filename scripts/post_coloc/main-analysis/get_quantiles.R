## Get quantiles of CLPP and CLPP-mod scores, so that
## we have a better idea of what a given cutoff means

## For the final analysis, we counted colocalizations as
## those that fell in the top 80% of tests, which agrees with
## our intuitive interpretation

require(qvalue)

coloc_results_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/full_coloc_results_qced.txt" 

# TODO: Fix path
results = read.table(coloc_results_file, sep="\t", header=TRUE)
results$gwas_p = 10^(-results$neg_log_gwas_pval)
results$eqtl_p = 10^(-results$neg_log_eqtl_pval)

sqtl = results[grepl("sQTL", results$eqtl_file),]
eqtl = results[grepl("eQTL", results$eqtl_file),]


# Quantiles for various cutoffs
allq = round(quantile(results$clpp_mod, seq(0,1,by=0.05)),4)
sq = round(quantile(sqtl$clpp_mod, seq(0,1,by=0.05)),4)
eq = round(quantile(eqtl$clpp_mod, seq(0,1,by=0.05)),4)

dat = data.frame(quant=seq(0,1,by=0.05), all=allq, sqtl=sq, eqtl=eq)
print(dat)

# Quantiles for various cutoffs
allq = round(quantile(results$clpp, seq(0,1,by=0.05)),4)
sq = round(quantile(sqtl$clpp, seq(0,1,by=0.05)),4)
eq = round(quantile(eqtl$clpp, seq(0,1,by=0.05)),4)

dat = data.frame(quant=seq(0,1,by=0.05), all=allq, sqtl=sq, eqtl=eq)
print(dat)



# Histogram of values
#hist(1-results$clpp)
hist(1-results$clpp_mod)


