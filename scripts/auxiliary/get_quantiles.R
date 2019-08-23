## Get quantiles of CLPP and CLPP-mod scores, so that
## we have a better idea of what a given cutoff means

require(qvalue)

results = read.csv("/users/mgloud/projects/insulin_resistance/output/full_coloc_results_qced.txt", sep=",")
results$gwas_p = 10^(-results$X.log_gwas_pval)
results$eqtl_p = 10^(-results$X.log_eqtl_pval)

# Quantiles for various cutoffs
quantile(results$clpp_mod, seq(0,1,by=0.05))
quantile(results[(results$X.log_gwas_pval > -log10(5e-8)) & (results$X.log_eqtl_pval > -log10(1e-5)),]$clpp_mod, seq(0,1,by=0.05))
quantile(results[(results$X.log_gwas_pval > -log10(1e-5)) & (results$X.log_eqtl_pval > -log10(1e-5)),]$clpp_mod, seq(0,1,by=0.05))
quantile(results[(results$X.log_gwas_pval > -log10(5e-8)) & (results$X.log_eqtl_pval > -log10(1e-6)),]$clpp_mod, seq(0,1,by=0.05))

filtered = results[(results$X.log_gwas_pval > -log10(5e-8)) & (results$X.log_eqtl_pval > -log10(1e-6)),]

# Histogram of values
hist(1-results$clpp)
hist(1-results$clpp_mod)

hist(1-filtered$clpp)
hist(1-filtered$clpp_mod)

qvalue(p=1-results$clpp)
qvalue(p=1-results$clpp_mod)
qvalue(p=1-filtered$clpp)
qvalue(p=1-filtered$clpp_mod)

