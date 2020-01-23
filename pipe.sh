
# Munge all files of interest
bash /users/mgloud/projects/insulin_resistance/scripts/format_gwas/munge.sh

# List SNPs to test
python /users/mgloud/projects/gwas-download/overlap/list_snps_to_test.py /users/mgloud/projects/insulin_resistance/scripts/config/ir.eqtls.only.overlap.config

# Run colocalization pipeline
python colocalization/test_IR_colocs.py

#
# (For eQTLs only)
#

# Make table with colocalization results.
# Get basic metrics on the number of runs and candidates;
# Also perform basic QC checks
Rscript /users/mgloud/projects/insulin_resistance/scripts/colocalization/render_output.R /users/mgloud/projects/insulin_resistance/scripts/colocalization/eqtls_only.post.config

# A few last QC sanity checks to make sure none of the colocalization tests have been dropped,
# as we've sometimes seen in the past. Consider merging these with the summarize_pre_test_statistics.Rmd file though
Rscript /users/mgloud/projects/insulin_resistance/scripts/colocalization/last_qc_check.R

# Bin coloc results into different categories of loci
# based on Ivan's method
Rscript /users/mgloud/projects/insulin_resistance/scripts/classify_results/classify_results.R /users/mgloud/projects/insulin_resistance/scripts/classify_results/ir_classification_2019-09-08.config

#
# (Repeat for sQTLs)
#

# Make table with colocalization results

Rscript /users/mgloud/projects/insulin_resistance/scripts/colocalization/render_output.R /users/mgloud/projects/insulin_resistance/scripts/colocalization/sqtls_only.post.config
Rscript /users/mgloud/projects/insulin_resistance/scripts/classify_results/classify_results.R /users/mgloud/projects/insulin_resistance/scripts/classify_results/ir_classification_sqtls_2019-09-08.config
# Bin coloc results into different categories of loci
# based on Ivan's method
