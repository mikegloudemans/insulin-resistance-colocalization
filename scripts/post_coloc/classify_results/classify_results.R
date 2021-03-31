# Author: Mike Gloudemans

#
# Implement Ivan's method for sorting out the test results
# into various categories
#
# 8/1/2019: Generalized it so that the three steps can now be
# run with arbitrary specifications, outlined in an accompanying
# config file.
#
# TODO: Provide documentation if we're going to distribute this
# publicly
#

require(dplyr)
require(rjson)

args = commandArgs(trailingOnly=TRUE)

# The config file for the whole post-coloc analysis, not just the classification step
config_file = args[1]

# TODO: The way I'm currently specifying the date of file creation is atrocious.
# Fix this to make it clean. 
full_config = fromJSON(file = config_file)

# TODO: Fix this to make the configs just be one file and one R object
# instead of having them be some weird multi-file thing
# 
# For now, the directory to the classification config file should
# be specified in the config for the entire pipeline run
class_config = fromJSON(file = full_config$classification_config_file)

# Load colocalization results matrix (output from the other
# QC and preparation script)
data = read.table(paste0(full_config$out_dir, "/clpp_results_", full_config$analysis_date, ".txt"), header=TRUE, sep="\t")

# Every SNP should be assigned to a locus!
stopifnot(sum(is.na(data$locus)) == 0)

cutoff_gwas_pval = as.numeric(class_config$pre_filters$max_gwas_pval)
cutoff_eqtl_pval = as.numeric(class_config$pre_filters$max_eqtl_pval)
cutoff_snp_count = as.numeric(class_config$pre_filters$min_snp_count)

if ("clpp_quantile_threshold" %in% names(class_config$pre_filters))
{
	clpp_threshold = quantile(as.numeric(data$clpp_mod), as.numeric(class_config$pre_filters$clpp_quantile_threshold)) 
} else
{
	clpp_threshold = as.numeric(class_config$pre_filters$clpp_threshold)
}

###################################
# Part 0: Pre-filtering
###################################

# Subset down to the ones that meet our initial screening criteria
# NOTE: This shouldn't really be necessary because if it passed the QC
# check then nothing should be different

# If they pass these p-values, they definitely pass the cut
eqtl_passing = (data$min_eqtl_pval < cutoff_eqtl_pval)    # eQTL significant
gwas_passing = (data$min_gwas_pval < cutoff_gwas_pval)    # GWAS significant

# But there are also a few exceptions...
# These GWAS get a boost
for (study in names(class_config$pre_filters$exceptions$gwas))
{
	extra_gwas = (grepl(study, data$base_gwas_file) &
		      (data$min_gwas_pval < as.numeric(class_config$pre_filters$exceptions$gwas[[study]]$max_gwas_pval)))
	gwas_passing = gwas_passing | extra_gwas
}
# These eQTLs get a boost
for (study in names(class_config$pre_filters$exceptions$eqtl))
{
	extra_eqtl = (grepl(study, data$base_gwas_file) &
		      (data$min_eqtl_pval < as.numeric(class_config$pre_filters$exceptions$eqtl[[study]]$max_eqtl_pval)))
	eqtl_passing = eqtl_passing | extra_eqtl
}

# Now filter down our matrix to those passing all pre-test thresholds.
sub = data[gwas_passing & eqtl_passing,]

# Remove those not having enough SNPs to pass filters
sub = sub[sub$n_snps >= cutoff_snp_count,]

sub$ref_snp = paste(sub$chr, sub$pos, sep="_")

loci_list = sort(unique(sub$locus))
step1_list = rep("", length(loci_list))
step2_list = rep("", length(loci_list))
step3_list = rep("", length(loci_list))
pancreas_list = rep("", length(loci_list))

###################################
# Part 1: General sorting
###################################


summary = sub %>% group_by(locus, ensembl) %>% summarize(colocs=sum(clpp_mod >= clpp_threshold))
coloc_counts = summary %>% group_by(locus) %>% summarize(coloc_genes = sum(colocs > 0), any_coloc_genes=sum(colocs > 0), candidates=length(ensembl))

candidates_equals = function(x, loci, coloc_matrix)
{
	return(loci[loci %in% (coloc_matrix[coloc_counts$candidates == x,]$locus)])
}
coloc_equals = function(x, loci, coloc_matrix)
{
	return(loci[loci %in% (coloc_matrix[coloc_counts$coloc_genes == x,]$locus)])
}
candidates_greater_than = function(x, loci, coloc_matrix)
{
	return(loci[loci %in% (coloc_matrix[coloc_counts$candidates > x,]$locus)])
}
coloc_greater_than = function(x, loci, coloc_matrix)
{
	return(loci[loci %in% (coloc_matrix[coloc_counts$coloc_genes > x,]$locus)])
}
for (type in names(class_config$step1_coloc_sorting))
{
	pass = coloc_counts$locus

	if ("candidates_equals" %in% names(class_config$step1_coloc_sorting[[type]][[1]]))
	{
		pass = candidates_equals(class_config$step1_coloc_sorting[[type]][[1]]["candidates_equals"], pass, coloc_counts)
	}
	if ("coloc_equals" %in% names(class_config$step1_coloc_sorting[[type]][[1]]))
	{
		pass = coloc_equals(class_config$step1_coloc_sorting[[type]][[1]]["coloc_equals"], pass, coloc_counts)
	}
	if ("candidates_greater_than" %in% names(class_config$step1_coloc_sorting[[type]][[1]]))
	{
		pass = candidates_greater_than(class_config$step1_coloc_sorting[[type]][[1]]["candidates_greater_than"], pass, coloc_counts)
	}
	if ("coloc_greater_than" %in% names(class_config$step1_coloc_sorting[[type]][[1]]))
	{
		pass = coloc_greater_than(class_config$step1_coloc_sorting[[type]][[1]]["coloc_greater_than"], pass, coloc_counts)
	}
	stopifnot(sum(step1_list[which(loci_list %in% pass)] != "") == 0)
	step1_list[which(loci_list %in% pass)] = type
}

# All loci should belong to a group at this point; if not, the groups are misspecified
# TODO: Allow for the possibility of "other" group for all not meeting these specs
stopifnot(sum(step1_list == "") == 0)
	
###################################
# Part 2: Tissue-specificity
###################################

# Figure out which tissues had colocs at each locus
tissue_coloc = sub %>% group_by(locus, eqtl_file) %>% summarize(has_coloc = as.numeric(sum(clpp_mod > clpp_threshold) > 0), has_no_coloc = as.numeric(sum(clpp_mod > clpp_threshold) == 0))

# Make sure all groups have been assigned to exactly one of these classes
stopifnot(sum(rowSums(tissue_coloc[,3:4]) == 0) == 0)

coloc = tissue_coloc[tissue_coloc$has_coloc == 1,] 
coloc_loci = unique(coloc$locus)
coloc_classes = sapply(coloc_loci, function(x)
       {

		tissues = coloc[coloc$locus == x,]$eqtl_file
       		# Test whether the colocalized tissues match some
       		# tissue group of interest

       		# See if there's a coloc in Pancreas
       		for (i in 1:length(class_config$step2_tissue_sorting))
		{
			type = class_config$step2_tissue_sorting[[i]]
			# It's a match if there are no tissues colocalized
			# that are not in this class
			if(sum(!(tissues %in% type$tissues)) == 0)
			{
				return(type$name)
			}

		}

		# If not, it still does have at least one coloc,
		# so doesn't belong in the "None" category
		return("Other")
       })

step2_list[which(loci_list %in% coloc_loci)] = coloc_classes

# Finally, the rest should have no colocs at all
# Make sure these loci actually have no colocs, then tag them
coloc_by_locus = tissue_coloc %>% group_by(locus) %>% summarize(total_coloc = sum(has_coloc) == 0)
stopifnot(sum(!(step2_list[which(loci_list %in% coloc_by_locus$locus[coloc_by_locus$total_coloc])] == "")) == 0)
step2_list[which(loci_list %in% coloc_by_locus$locus[coloc_by_locus$total_coloc])] = "None"

pancreas_coloc = tissue_coloc[grepl("Pancreas", tissue_coloc$eqtl_file) & (tissue_coloc$has_coloc == 1),]
pancreas_list[rep(TRUE, length(pancreas_list))] = "none"
pancreas_list[which(loci_list %in% unique(pancreas_coloc$locus))] = "coloc"


# Just make sure every site's been assigned now
stopifnot(sum(step2_list == "") == 0)


###################################
# Part 3: Which GWAS matter most?
###################################

top_colocs = sub %>% group_by(locus, base_gwas_file) %>% summarize(best = max(clpp_mod))

# Tag colocalized loci by priority level
gwas_cumul_loci = c()
for (i in 1:length(class_config$step3_gwas_sorting))
{
	gwas_loci = unique(top_colocs[(top_colocs$base_gwas_file %in% class_config$step3_gwas_sorting[[i]][["traits"]]) & (top_colocs$best > clpp_threshold),]$locus)
	gwas_loci = gwas_loci[!(gwas_loci %in% gwas_cumul_loci)]
	step3_list[which(loci_list %in% gwas_loci)] = class_config$step3_gwas_sorting[[i]][["name"]]

	gwas_cumul_loci = c(gwas_cumul_loci, gwas_loci)
}

# The remainder of colocalizing loci will be tagged as "other".
# This only matters if we have traits that are in none of the tiers.
loci_other = top_colocs[top_colocs$best > clpp_threshold,]$locus
loci_other = loci_other[!(loci_other %in% gwas_cumul_loci)]
step3_list[which(loci_list %in% loci_other)] = "Other"

# Make sure all the remaining loci don't actually have any colocalization whatsoever
# (as computed for step 2 of prioritization)
stopifnot(sum(!(step3_list[which(loci_list %in% coloc_by_locus$locus[coloc_by_locus$total_coloc])] == "")) == 0)
step3_list[which(loci_list %in% coloc_by_locus$locus[coloc_by_locus$total_coloc])] = "None"
# Now just make sure all have been tagged
stopifnot(sum(step3_list == "") == 0)


#############################################
# Get names of colocalized genes. We want
# - all candidate eQTL genes
# - all candidate sQTL genes
# - all candidate sQTL features
# - all candidate genes (either)
# - all colocalized eQTL genes
# - all colocalized sQTL genes
# - all colocalized sQTL features
# - all colocalized genes (either)
#############################################
coloc_genes = sub %>% group_by(locus) %>% summarize(candidate_genes = paste(unique(ensembl), collapse=";"), 
				      candidate_eqtl_genes = paste(unique(ensembl[grepl("eQTL", eqtl_file)]), collapse=";"), 
				      candidate_sqtl_genes = paste(unique(ensembl[grepl("sQTL", eqtl_file)]), collapse=";"), 
				      candidate_sqtl_features = paste(unique(feature[grepl("sQTL", eqtl_file)]), collapse=";"),
				      coloc_genes = paste(unique(ensembl[(clpp_mod >= clpp_threshold)]), collapse=";"), 
				      coloc_eqtl_genes = paste(unique(ensembl[(clpp_mod >= clpp_threshold)& grepl("eQTL", eqtl_file)]), collapse=";"), 
				      coloc_sqtl_genes = paste(unique(ensembl[(clpp_mod >= clpp_threshold)& grepl("sQTL", eqtl_file)]), collapse=";"), 
				      coloc_sqtl_features = paste(unique(feature[(clpp_mod >= clpp_threshold) & grepl("sQTL", eqtl_file)]), collapse=";"),
				      coloc_both = paste(unique(ensembl[(clpp_mod >= clpp_threshold)& grepl("eQTL", eqtl_file)])[unique(ensembl[(clpp_mod >= clpp_threshold)& grepl("eQTL", eqtl_file)]) %in% unique(ensembl[(clpp_mod >= clpp_threshold)& grepl("sQTL", eqtl_file)])], collapse=";"),
				      coloc_traits = paste(unique(gwas_short[(clpp_mod >= clpp_threshold)]), collapse=","),
				      all_lead_gwas_snps=paste(unique(rsid[(clpp_mod >= clpp_threshold)]),collapse=","),
				      all_lead_gwas_coords=paste(unique(ref_snp[(clpp_mod >= clpp_threshold)]),collapse=",")
				      )


####### Put all the results together ########

classes = data.frame(list(locus=loci_list, step1=step1_list, step2=step2_list, step3=step3_list, pancreas=pancreas_list))

classes = merge(classes, coloc_genes)

write.table(classes, file=paste0(full_config$out_dir, "/coloc_classification_", full_config$analysis_date, ".txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


results_summaries = list()
results_summaries$step1 = classes %>% group_by(step1) %>% summarize(length(locus))
results_summaries$step2 = classes %>% group_by(step2) %>% summarize(length(locus))
results_summaries$step3 = classes %>% group_by(step3) %>% summarize(length(locus))
results_summaries$pancreas = classes %>% group_by(pancreas) %>% summarize(length(locus))


summary_file = paste0(full_config$out_dir, "/coloc_class_summary_", full_config$analysis_date, ".txt")

# Remove this file if it already exists
suppressWarnings(file.remove(file=summary_file))
suppressWarnings(lapply(results_summaries, function(x) {write.table( data.frame(x), summary_file, append= T, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE); write("\n\n",file=summary_file, append=TRUE)}))

# Add classifications to the original data frame and rewrite it to a file.

data_extended = full_join(data, classes)
write.table(data_extended, file=paste0(full_config$out_dir, "/clpp_results_categorized_", full_config$analysis_date, ".txt"), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

