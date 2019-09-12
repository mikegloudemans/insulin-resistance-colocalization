
# Munge all files of interest
bash /users/mgloud/projects/insulin_resistance/scripts/format_gwas/munge.sh

# List SNPs to test
python /users/mgloud/projects/gwas-download/overlap/list_snps_to_test.py /users/mgloud/projects/insulin_resistance/scripts/config/ir.overlap.config

# Run colocalization pipeline
python colocalization/test_IR_colocs.py

#
# (For eQTLs only)
#

# Make table with colocalization results

Rscript /users/mgloud/projects/insulin_resistance/scripts/colocalization/render_output.R /users/mgloud/projects/insulin_resistance/scripts/colocalization/eqtls_only.post.config
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
