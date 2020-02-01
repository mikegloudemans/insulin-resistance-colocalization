# Full pipeline for running all colocalization and post-coloc categorization
# Author: Mike Gloudemans
# Updated 1/23/2020

# Munge all files of interest
bash scripts/format_gwas/munge.sh

# List SNPs to test
python bin/gwas-download/overlap/list_snps_to_test.py scripts/pre_coloc/overlap_config/ir.overlap.config

# Use this line instead to run for only eQTLs and no other QTL types
#python bin/gwas-download/overlap/list_snps_to_test.py scripts/pre_coloc/overlap_config/ir.eqtls.only.overlap.config

# Run colocalization pipeline
python scripts/colocalization/test_IR_colocs.py

#
# (For eQTLs only)
#

# Make table with colocalization results.
# Get basic metrics on the number of runs and candidates;
# Also perform basic QC checks
Rscript scripts/post_coloc/summarize_pre_test_statistics.R scripts/post_coloc/eqtls_only.post.config
Rscript scripts/post_coloc/aggregate_coloc_results.R scripts/post_coloc/eqtls_only.post.config

# A few last QC sanity checks to make sure none of the colocalization tests have been dropped,
# as we've sometimes seen in the past. Consider merging these with the summarize_pre_test_statistics.Rmd file though
#
# NOTE: This is a "soft QC" step; it will not automatically throw an error if the results are
# improper. The user should visually inspect the printed output to verify that not many tests
# were inexplicably dropped, no tests passed through p-value filters, every GWAS is represented, etc.
Rscript scripts/post_coloc/last_qc_check.R

# Bin coloc results into different categories of loci
# based on Ivan's method
Rscript scripts/classify_results/classify_results.R scripts/post_coloc/ir_classification_2019-09-08.config

#
# (Repeat for sQTLs)
#

# Make table with colocalization results

Rscript scripts/post_coloc/render_output.R scripts/post_coloc/sqtls_only.post.config
Rscript scripts/post_coloc/classify_results/classify_results.R scripts/post_coloc/classify_results/ir_classification_sqtls_2019-09-08.config
# Bin coloc results into different categories of loci
# based on Ivan's method

#########################################################################
### Additional analysis and checks that were used but not explicitly
### necessary to generate the main results
#########################################################################

# Estimating a typical CLPP or CLPP-mod score at various quantiles
scripts/auxiliary/get_quantiles.R



