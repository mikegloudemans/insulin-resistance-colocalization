# Full pipeline for running all colocalization and post-coloc categorization
# Author: Mike Gloudemans

####################################################
# Make the directories that we'll need later
# Anything deeper will be created by the individual 
# scripts at runtime
####################################################

mkdir tmp

mkdir -p output/colocalization
mkdir -p output/post_coloc
mkdir -p output/test_snps

####################################################
# Pre-processing steps to format the relevant data 
####################################################

# Munge all files of interest
bash scripts/format_gwas/munge.sh
# ^ Admittedly this is missing the config file for the newest T2D files, which were done with
# a separate script in the munge module

# List SNPs to test
#python bin/gwas-download/overlap/list_snps_to_test.py scripts/pre_coloc/overlap_config/ir.overlap.config
# ^ obsolete now that t2d.extended exists
python bin/gwas-download/overlap/list_snps_to_test.py scripts/pre_coloc/overlap_config/ir.overlap.t2d.extended.config

# Use this line instead to run for only eQTLs and no other QTL types
#python bin/gwas-download/overlap/list_snps_to_test.py scripts/pre_coloc/overlap_config/ir.eqtls.only.overlap.config

# exclude edQTLs for now
grep -v aggro output/test_snps/ir-v8/ir-v8-extended_all-gwas_gtex-ir_source-pval1e-05_lookup-pval1e-05_source-window500000_lookup-window10000_coloc-tests.txt > output/test_snps/ir-v8/ir-v8-extended-coloc-tests-no-edQTL.txt

# Run colocalization pipeline for all QTL types
python bin/colocalization_pipeline/dispatch.py scripts/colocalization/ir-t2d-extended.config

#
# A general post-processing loop that can be applied for any QTL type
# given the correct config file
#

eqtl_config_file="scripts/post_coloc/config/eqtls_only.post.config"
sqtl_config_file="scripts/post_coloc/config/sqtls_only.post.config"
sqtl_or_eqtl_config_file="scripts/post_coloc/config/eqtls_and_sqtls.post.config"

# Just add more files to the list if we want to repeat post-processing for other QTL types. The steps are the
# same, as long as the appropriate config file is used.
for post_config_file in $eqtl_config_file $sqtl_config_file $sqtl_or_eqtl_config_file; do
	echo $post_config_file

	# Make table with colocalization results.
	# Get basic metrics on the number of runs and candidates;
	# Also perform basic QC checks
	Rscript scripts/post_coloc/summarize_pre_test_statistics.R $post_config_file
	Rscript scripts/post_coloc/aggregate_coloc_results.R $post_config_file

	# A few last QC sanity checks to make sure none of the colocalization tests have been dropped,
	# as we've sometimes seen in the past. Consider merging these with the summarize_pre_test_statistics.Rmd file though
	#
	# NOTE: This is a "soft QC" step; it will not automatically throw an error if the results are
	# improper. The user should visually inspect the printed output to verify that not many tests
	# were inexplicably dropped, no tests passed through p-value filters, every GWAS is represented, etc.
	Rscript scripts/post_coloc/last_qc_check.R $post_config_file

	# Bin coloc results into different categories of loci
	# based on Ivan's method
	Rscript scripts/post_coloc/classify_results/classify_results.R $post_config_file

	Rscript scripts/post_coloc/get_variant_veps.R $post_config_file 
done

#########################################################################
### Additional analysis and checks that were used but not explicitly
### necessary to generate the main results
#########################################################################

# Estimating a typical CLPP or CLPP-mod score at various quantiles
scripts/auxiliary/get_quantiles.R


#############################
# In progress
#############################

