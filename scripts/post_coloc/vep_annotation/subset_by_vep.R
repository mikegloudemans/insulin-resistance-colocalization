require(dplyr)

veps = read.table("output/post_coloc/2019-12-20/refiltered/eqtls_only/vep_summary_2019-12-20.txt", header=FALSE, skip=1)
melt_vep = read.table("output/post_coloc/2019-12-20/refiltered/eqtls_only/vep_summary_melted_2019-12-20.txt", header=TRUE)
categories = read.table("output/post_coloc/2019-12-20/refiltered/eqtls_only/clpp_results_categorized_2019-12-20.txt", header=TRUE, fill=TRUE, sep="\t")
names(veps) =c("locus", "veps")


test_set_for_veps = function(strong_coloc)
{
	full_coloc = merge(strong_coloc, veps, by="locus")
	length(unique(full_coloc$locus))

	vep_types = unique(melt_vep$vep)

	# How many of our loci of interest have each VEP consequence?
	for (consequence in vep_types)
	{
		if(is.na(consequence)){next}

		print(consequence)
		loci_with_effect = full_coloc %>% group_by(locus) %>% summarize(has_consequence = sum(sapply(veps, function(x) {grepl(consequence, x)})) > 0)
		print(sum(loci_with_effect$has_consequence))
	}

	# How many of them have any of the obvious high-risk consequences?
	high_risk_summary = full_coloc %>% group_by(locus) %>% summarize(has_consequence = sum(sapply(veps, function(x) {grepl("NON_SYNONYMOUS", x)})
											       + sapply(veps, function(x) {grepl("CANONICAL_SPLICE", x)})
											       + sapply(veps, function(x) {grepl("STOP_GAINED", x)}))>0)
	print("Number of loci with a high-risk consequence:")
	print(sum(high_risk_summary$has_consequence))
	print("Number of loci without a high-risk consequence:")
	print(sum(!high_risk_summary$has_consequence))

	alt_scaffold = full_coloc %>% group_by(locus) %>% summarize(has_alt = sum(sapply(veps, function(x) {grepl("alt", x)})
											  + sapply(veps, function(x) {grepl("alt", x)}))>0)
	print("Number of loci in LD with an alternate genome scaffold:")
	print(sum(alt_scaffold$has_alt))
}

print("All loci:")
print("")
print("TOTAL:")
print(length(unique(categories$locus)))
print("")
test_set_for_veps(categories)

print("")
print("All strong coloc loci:")
strong_coloc = categories[categories$step1 %in% c("loci_0.1", "loci_2", "loci_3"),]
print("")
print("TOTAL:")
print(length(unique(strong_coloc$locus)))
print("")
test_set_for_veps(strong_coloc)

print("")
print("All strong coloc loci with only one coloc gene:")
strong_coloc = categories[categories$step1 %in% c("loci_0.1", "loci_2"),]
print("")
print("TOTAL:")
print(length(unique(strong_coloc$locus)))
print("")
test_set_for_veps(strong_coloc)
