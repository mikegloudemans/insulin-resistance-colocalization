# Compare list of tests to run and list of tests actually run.
# There shouldn't be many differences.
# For the differences that exist, is there a reasonable explanation?

pre_run_coloc = read.table("/users/mgloud/projects/insulin_resistance/output/test_snps/ir-v8_all-gwas_gtex-ir_gwas-pval1e-05_eqtl-pval1e-05_gwas-window500000_eqtl-window0_coloc-tests.txt", header = TRUE)
post_run_coloc = read.table("/users/mgloud/projects/insulin_resistance/output/post_coloc/eqtls_only/full_coloc_results_qced.txt", header =TRUE)

pre_run_coloc$test_id = paste(pre_run_coloc$chr, pre_run_coloc$snp_pos, pre_run_coloc$eqtl_file, pre_run_coloc$gwas_file, pre_run_coloc$feature, sep="_")
post_run_coloc$test_id = paste(post_run_coloc$ref_snp, post_run_coloc$eqtl_file, post_run_coloc$base_gwas_file, post_run_coloc$feature, sep="_")

# Make this more general at some point
pre_run_coloc = pre_run_coloc[!grepl("2hr", pre_run_coloc$gwas_file),]
pre_run_coloc = pre_run_coloc[(grepl("ISI_Model_", pre_run_coloc$gwas_file) | grepl("MI_adjBMI_", pre_run_coloc$gwas_file) | (pre_run_coloc$gwas_pvalue < 5e-8)),]

sum(pre_run_coloc$test_id %in% post_run_coloc$test_id)
sum(!(pre_run_coloc$test_id %in% post_run_coloc$test_id))

rejected_tests = pre_run_coloc[!(pre_run_coloc$test_id %in% post_run_coloc$test_id),]


# A few other QC checks

# We shouldn't see many sites with < 50 or so SNPs. We shouldn't really see any
# with < 10 SNPs; if we do, we need to toss these sites at an earlier step in the
# process.
# TODO: Add a min_snps parameter to the pre-run config file
print(range(post_run_coloc$n_snps))

# The violators are all FastGlu and FastInsu SNPs, which isn't surprising.
# This is a good sanity check though.
print(post_run_coloc[post_run_coloc$n_snps < 50,]$base_gwas_file)

