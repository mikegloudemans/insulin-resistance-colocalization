# Compare list of tests to run and list of tests actually run.
# There shouldn't be many differences.
# For the differences that exist, is there a reasonable explanation?

suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(rjson)))

config_file = commandArgs(trailingOnly=TRUE)[1]

# Load pre-specified config file
config = fromJSON(file=config_file)

pre_run_coloc = read.table(config$coloc_test_list, header=TRUE)
post_run_coloc = read.table(paste(config$out_dir, "full_coloc_results_qced.txt", sep="/"), header=TRUE, sep="\t")

pre_run_coloc$test_id = paste(pre_run_coloc$chr, pre_run_coloc$snp_pos, pre_run_coloc$lookup_file, pre_run_coloc$source_file, pre_run_coloc$lookup_trait, sep="_")
post_run_coloc$test_id = paste(post_run_coloc$ref_snp, post_run_coloc$eqtl_file, post_run_coloc$base_gwas_file, post_run_coloc$feature, sep="_")

# TODO: Make this more general at some point
# Right now it's not reflecting the filters we applied at the last step
pre_run_coloc = pre_run_coloc[!grepl("2hr", pre_run_coloc$source_file),]
pre_run_coloc = pre_run_coloc[(grepl("ISI_Model_", pre_run_coloc$source_file) | grepl("MI_adjBMI_", pre_run_coloc$source_file) | (pre_run_coloc$source_pvalue < 5e-8)),]

print("Number of candidate loci successfully tested:")
sum(pre_run_coloc$test_id %in% post_run_coloc$test_id)
print("")
print("Number of candidate loci not appearing in final results:")
sum(!(pre_run_coloc$test_id %in% post_run_coloc$test_id))
print("")

rejected_tests = pre_run_coloc[!(pre_run_coloc$test_id %in% post_run_coloc$test_id),]

# Make sure the range of GWAS p-values for each type of test is reasonable

# Visual inspection:
# Each GWAS should have a reasonable number of tested SNPs.
# The range of SNP pvalues for each GWAS and eQTL study should match the specified cutoffs
print("Number of tested SNPs and p-value range for each GWAS:")
gwas_stats = post_run_coloc %>% group_by(base_gwas_file) %>% summarize(total_snps=length(unique(ref_snp)), gwas_range=paste(round(range(neg_log_gwas_pval),2), collapse="..."), eqtl_range = paste(round(range(neg_log_eqtl_pval),2), collapse="..."))
print(gwas_stats)
print("")

# Same as above, just grouped by eQTLs instead of GWAS
print("Number of tested SNPs and p-value range for each QTL study:")
eqtl_stats = post_run_coloc %>% group_by(eqtl_file) %>% summarize(total_snps=length(unique(ref_snp)), gwas_range=paste(round(range(neg_log_gwas_pval),2), collapse="..."), eqtl_range = paste(round(range(neg_log_eqtl_pval),2), collapse="..."))
print(eqtl_stats)
print("")

# We shouldn't see many sites with < 50 or so SNPs. We shouldn't really see any
# with < 10 SNPs; if we do, we need to toss these sites at an earlier step in the
# process.
# TODO: Add a min_snps parameter to the pre-run config file
print("Range of number of SNPs per locus:")
print(range(post_run_coloc$n_snps))
print("")

# The violators are all FastGlu and FastInsu SNPs, which isn't surprising.
# This is a good sanity check though.
print("Loci with low number of SNPs:")
print(table(post_run_coloc[post_run_coloc$n_snps < 50,]$base_gwas_file))
print("")

