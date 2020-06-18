# Some checks that should be run at the end of summarize_pre_test_statistics.R if
# we're trying to entire "fill the matrix" i.e. run every tissue x gene x locus x GWAS


if (FALSE)

{

	#########################################################
	### QC stuff: troubleshooting failed tests
	#########################################################

	### How many SNPs were dropped during the colocalization analysis?
	dropped = testable_snps[!(testable_snps$ref_snp %in% unique(results$ref_snp)),]
	print("Tests dropped:")
	dim(dropped)[1]
	print("SNP-gene pairs dropped:")
	sum(!duplicated(dropped[c("ref_snp", "feature")]))
	print("SNPs dropped:")
	length(unique(dropped$ref_snp))

	### Which SNPs were tested even though they don't appear to pass the required threshold in any tissue?

	# NOTE: I think this section only matters if we're running the mode where we fill the matrix, allowing all
	# SNPs to be tested if they end up significant in any trait whatsoever. But the above filters may be problematic
	# in such a case...the bottom line is that for the final generalizable version, I should choose one as the
	# default and then add the other option (probably "matrix-filling") as a user-definable parameter

	# Why? (Since we apply a filter at the beginning, this suggests that
	# they've been dropped because their lead variant was not measured in the eQTL study.)

	best_pvals = results %>% group_by(ref_snp, feature) %>% summarize(best_gwas_pval = max(X.log_gwas_pval), best_eqtl_pval = max(X.log_eqtl_pval))

	best_pvals$best_gwas = sapply(1:dim(best_pvals)[1], function(i)
	       {
			gene = best_pvals$feature[i]
			pval = best_pvals$best_gwas_pval[i]
			snp = best_pvals$ref_snp[i]

			best_pval = paste(unique(results[results$ref_snp == snp & abs(results$X.log_gwas_pval - pval) <= 0.01 & results$feature == gene,]$base_gwas_file), collapse="-")


			return(best_pval)
	       })

	best_pvals$best_eqtl = sapply(1:dim(best_pvals)[1], function(i)
	       {
			gene = best_pvals$feature[i]
			pval = best_pvals$best_eqtl_pval[i]
			snp = best_pvals$ref_snp[i]

			best_pval = paste(unique(results[results$ref_snp == snp & abs(results$X.log_eqtl_pval - pval) <= 0.01 & results$feature == gene,]$eqtl_file), collapse="-")
			return(best_pval)
	       })

	# For later, we'll also want to know all traits in which the GWAS and the eQTL p-values were significant.
	results$all_sig_gwas = sapply(1:dim(results)[1], function(i)
		{
			feature = results$feature[i]
			ref_snp = results$ref_snp[i]
			sub = results[(results$feature == feature) & (results$ref_snp == ref_snp),]

			gwas_sig = sub[apply(sub[c("gwas_pval", "base_gwas_file")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]
			
			return(paste(unique(gwas_sig$base_gwas_file), collapse=";"))
		}
	)

	results$all_sig_eqtl = sapply(1:dim(results)[1], function(i)
		{
			feature = results$feature[i]
			ref_snp = results$ref_snp[i]
			sub = results[(results$feature == feature) & (results$ref_snp == ref_snp),]
			eqtl_sig = sub[apply(sub[c("eqtl_pval", "eqtl_file")], 1, FUN=pval_passing, threshold = config$eqtl_pval_threshold),]
			return(paste(unique(eqtl_sig$eqtl_file), collapse=";"))
		}
	)

	### Get the list of SNPs that shouldn't have actually been tested 

	# For these, the lead eQTL or GWAS SNP
	# was dropped during the intersection phase.

	best_pvals$best_eqtl_pval_true = 10^(-best_pvals$best_eqtl_pval)
	best_pvals$best_gwas_pval_true = 10^(-best_pvals$best_gwas_pval)

	insignificant_eqtl = best_pvals[!apply(best_pvals[c("best_eqtl_pval_true", "best_eqtl")], 1, FUN=pval_passing, threshold = config$eqtl_pval_threshold),]
	insignificant_gwas = best_pvals[!apply(best_pvals[c("best_gwas_pval_true", "best_gwas")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]

	print(dim(insignificant_eqtl)[1])
	print(insignificant_eqtl)
	print(dim(insignificant_gwas)[1])
	print(insignificant_gwas)

	# Note: If we want to follow up on a particular gene to see why it was dropped, we
	# can do it like so:
	# one_insig = results[(results$feature == "ENSG00000002919.10") & (results$ref_snp == "17_46292923"),]

	### Filter our results table to exclude the SNPs that are no longer significant.

	print("Results table dimension:")
	dim(results)[1]
	results = results[!(paste(results$ref_snp, results$feature, sep="_") %in% paste(insignificant_eqtl$ref_snp, insignificant_eqtl$feature, sep="_")),]
	results = results[!(paste(results$ref_snp, results$feature, sep="_") %in% paste(insignificant_gwas$ref_snp, insignificant_gwas$feature, sep="_")),]
	print("Results table dimension after filtering insignificant SNPs:")
	dim(results)[1]

	# Now that we have the final list of colocalization results,
	# we can assign a unique locus number to each locus, grouping
	# nearby SNPs.
	results$locus = group_to_loci(results$ref_snp)

	### Which SNPs were not tested for all trait-tissue combos? 

	# Why? (Maybe due to errors while running the pipeline -- see if these errors
	# are excusable.)

	tests_per_pair = results %>% group_by(ref_snp, feature) %>% summarize(tests = length(feature))
	missing_pairs = tests_per_pair[tests_per_pair$tests < (length(unique(results$eqtl_file)) * length(unique(results$base_gwas_file))),]
	print("Total number of snp-gene pairs tested:")
	dim(tests_per_pair)[1]
	print("Total number of snp-gene pairs missing at least one trait-tissue combo:")
	dim(missing_pairs)[1]
	all_missing_tests = results[paste(results$ref_snp, results$feature, sep="_") %in% paste(missing_pairs$ref_snp, missing_pairs$feature, sep="_"),]

	#I manually inspected the tests that were dropped. Based on the number of tests missed
	#for each locus, it's clear that dropped SNPs either occurred because the gene was not
	#tested for eQTLs in every trait, or because the GWAS summary stats at that locus had no overlap with the
	#eQTL summary statistics for at least one of the GWAS analyses.

	# For a bit more detail...

	### Why are tests missing for some SNP-gene pairs?

	# For snp-gene pairs that have missing tests, is it because there was an error thrown in the pipeline,
	# or because the variant was intentionally skipped due to non-overlap or being on the edge of the range?

	missing_pairs$error = paste(missing_pairs$ref_snp, missing_pairs$feature, sep="_") %in% unique(paste(errors$snp.chrom, errors$snp.pos, errors$restrict_gene, sep="_"))
	missing_pairs$skip = paste(missing_pairs$ref_snp, missing_pairs$feature, sep="_") %in% unique(paste(skips$snp.chrom, skips$snp.pos, skips$feature, sep="_"))
	sum(missing_pairs$error)
	sum(missing_pairs$skip)
	# Reasons why some combinations of trait-tissue at SNP-gene combos in our test set were skipped
	table(skips$error)

	# Based on the QC checks here, I'm confident that we're not missing any tests for unexplainable reasons.

	# After performing all of these checks, we're ready to go on to the main colocalization analysis.

	if (FALSE)
	{
		missing_pairs$reason_missing = sapply(1:dim(missing_pairs)[1], function(i)
		       {
				snp = missing_pairs$ref_snp[i]
				gene = missing_pairs$feature[i]
				sub = results[(results$ref_snp == snp) & (results$feature == gene),]
				
				code = 0
				# TODO: If we restore this code for "full matrix" form,
				# then we should substitute the numbers for something in config file
				if(length(unique(sub$eqtl_file)) < 4)
				{
					code = code + 1
				}
				if(length(unique(sub$gwas_trait)) < 9)
				{
					code = code + 2
				}
				return(code)

		       })

		write.table(missing_pairs[c("ref_snp", "feature", "reason_missing")], file=paste(config$out_dir, "all_missing_pairs_explained.txt", sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	}
}
