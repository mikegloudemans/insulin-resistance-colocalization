
# Munge all files of interest
bash /users/mgloud/projects/insulin_resistance/scripts/format_gwas/munge.sh

# List SNPs to test
python /users/mgloud/projects/gwas-download/overlap/list_snps_to_test.py /users/mgloud/projects/insulin_resistance/scripts/config/ir.overlap.config

# Run colocalization pipeline
python colocalization/test_IR_colocs.py

# Make table with colocalization results

# Bin coloc results into different categories of loci
# based on Ivan's method
