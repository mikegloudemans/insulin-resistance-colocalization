# Implement Ivan's method for sorting out the test results
# into various categories

require(dplyr)

# Load colocalization results matrix (output from the other
# QC and preparation script)
data = read.table("/users/mgloud/projects/insulin_resistance/scripts/classify_results/data/clpp_results_20190401.txt", header=TRUE, sep="\t")

# TODO: Make this general and load from a JSON config file

# In the general file, we'll have a field called "exceptions"
# for this sort of thing, I think

cutoff_gwas_pval = 5e-8
cutoff_eqtl_pval = 1e-5
cutoff_snp_count = 20

strong_clpp_threshold = 0.4
weak_clpp_threshold = 0.25

###################################
# Part 0: Pre-filtering
###################################

# Subset down to the ones that meet our initial screening criteria

# Remove all tests not meeting p-value cutoffs (unless specifically
# allowed to not meet those cutoffs)
sub = data[((data$min_eqtl_pval < cutoff_eqtl_pval) &		# eQTL significant
	    (data$min_gwas_pval < cutoff_gwas_pval)) 		# GWAS significant
		| ((((grepl("MI_adjBMI", data$gwas_trait) 	# Unless it's an exception trait
		      & data$min_gwas_pval < 1e-5) 
		| (grepl("ISI", data$gwas_trait) 
		      & (data$min_gwas_pval < 1e-5))))),]


# Remove those not having enough SNPs to pass filters
sub = sub[sub$n_snps >= cutoff_snp_count,]

loci_list = sort(unique(sub$locus))
step1_list = rep("", length(loci_list))
step2_list = rep("", length(loci_list))
step3_list = rep("", length(loci_list))

###################################
# Part 1: General sorting
###################################

summary = sub %>% group_by(locus, ensembl, hgnc) %>% summarize(colocs=sum(clpp_mod >= strong_clpp_threshold), weak_colocs=sum((clpp_mod >= weak_clpp_threshold) & (clpp_mod < strong_clpp_threshold)))
coloc_counts = summary %>% group_by(locus) %>% summarize(strong_coloc_genes = sum(colocs > 0), weak_coloc_genes = sum((weak_colocs > 0) & (colocs == 0)), any_coloc_genes=sum((colocs > 0) | (weak_colocs > 0)), candidates=length(ensembl))

# Category 0: Just one candidate, but at least it's colocalized
loci_0.1 = coloc_counts[(coloc_counts$candidates == 1) & (coloc_counts$strong_coloc_genes == 1),]$locus
loci_0.2 = coloc_counts[(coloc_counts$candidates == 1) & (coloc_counts$weak_coloc_genes == 1),]$locus

# Category 1: One candidate gene, and it doesn't even colocalize
loci_1 = coloc_counts[(coloc_counts$candidates == 1) & (coloc_counts$weak_coloc_genes == 0) & (coloc_counts$strong_coloc_genes == 0),]$locus

# Category 2: Multiple candidates, now narrowed down to one clear signal
loci_2 = coloc_counts[(coloc_counts$candidates > 1) & (coloc_counts$strong_coloc_genes == 1),]$locus

# Category 3: Multiple candidates, but more than one strong coloc signal
loci_3 = coloc_counts[(coloc_counts$candidates > 1) & (coloc_counts$strong_coloc_genes > 1),]$locus

# Category 4: Multiple candidates, any number of weak coloc signals
loci_4.1 = coloc_counts[(coloc_counts$candidates > 1) & (coloc_counts$strong_coloc_genes == 0) & (coloc_counts$weak_coloc_genes == 1),]$locus
loci_4.2 = coloc_counts[(coloc_counts$candidates > 1) & (coloc_counts$strong_coloc_genes == 0) & (coloc_counts$weak_coloc_genes > 1),]$locus

# Category 5: Multiple candidates, and none of them colocalize even weakly
loci_5 = coloc_counts[(coloc_counts$candidates > 1) & (coloc_counts$any_coloc_genes == 0),]$locus

# QC check: each locus should appear once and only once.
#
# NOTE: A few loci fell out already due to our thresholding steps,
# which is okay.
#
all_loci = c(loci_0.1, loci_0.2, loci_1, loci_2, loci_3, loci_4.1, loci_4.2, loci_5)
stopifnot(dim(coloc_counts)[1] == length(all_loci))
stopifnot(length(table(table(all_loci))) == 1)

step1_list[loci_list %in% loci_0.1] = "0.1"
step1_list[loci_list %in% loci_0.2] = "0.2"
step1_list[loci_list %in% loci_1] = "1"
step1_list[loci_list %in% loci_2] = "2"
step1_list[loci_list %in% loci_3] = "3"
step1_list[loci_list %in% loci_4.1] = "4.1"
step1_list[loci_list %in% loci_4.2] = "4.2"
step1_list[loci_list %in% loci_5] = "5"

###################################
# Part 2: Tissue-specificity
###################################

# Figure out which tissues had strong, weak, no colocs at each locus
tissue_coloc = sub %>% group_by(locus, eqtl_file) %>% summarize(has_strong_coloc = as.numeric(sum(clpp_mod > strong_clpp_threshold) > 0), has_weak_only = as.numeric((sum(clpp_mod > strong_clpp_threshold) == 0) & (sum(clpp_mod > weak_clpp_threshold) > 0)), has_no_coloc = as.numeric(sum(clpp_mod > weak_clpp_threshold) == 0))

# Make sure all groups have been assigned to exactly one of these classes
stopifnot(sum(rowSums(tissue_coloc[,3:5]) == 0) == 0)

# Get loci with no colocalization whatsoever, in any tissue
coloc_by_locus = tissue_coloc %>% group_by(locus) %>% summarize(total_coloc = sum(has_strong_coloc + has_weak_only) == 0)
loci_non_coloc = coloc_by_locus$locus[coloc_by_locus$total_coloc]
# Make sure this group includes the ones from groups 1 and 5 above
stopifnot(sum(loci_non_coloc %in% c(loci_1, loci_5) == 0) == 0)
stopifnot(sum(c(loci_1, loci_5) %in% loci_non_coloc == 0) == 0)

# First, let's just get the loci that have coloc in multiple tissues
strong_tissue_coloc = tissue_coloc[tissue_coloc$has_strong_coloc == 1,]
weak_tissue_coloc = tissue_coloc[tissue_coloc$has_weak_only == 1,]

strong_tissue_counts = strong_tissue_coloc %>% group_by(locus) %>% summarize(how_many = length(has_strong_coloc))
weak_tissue_counts = weak_tissue_coloc %>% group_by(locus) %>% summarize(how_many = length(has_weak_only))

loci_shared = c(strong_tissue_counts$locus[strong_tissue_counts$how_many > 1], weak_tissue_counts$locus[(weak_tissue_counts$how_many > 1) & !(weak_tissue_counts$locus %in% strong_tissue_counts$locus)])

adps_strong = strong_tissue_coloc$locus[strong_tissue_coloc$eqtl_file == "Adipose_Subcutaneous"]
adps_weak = weak_tissue_coloc$locus[weak_tissue_coloc$eqtl_file == "Adipose_Subcutaneous"]
loci_adps = c(adps_strong[!adps_strong %in% loci_shared], adps_weak[!adps_weak %in% c(loci_shared, strong_tissue_counts$locus)])

adpv_strong = strong_tissue_coloc$locus[strong_tissue_coloc$eqtl_file == "Adipose_Visceral_Omentum"]
adpv_weak = weak_tissue_coloc$locus[weak_tissue_coloc$eqtl_file == "Adipose_Visceral_Omentum"]
loci_adpv = c(adpv_strong[!adpv_strong %in% loci_shared], adpv_weak[!adpv_weak %in% c(loci_shared, strong_tissue_counts$locus)])

musk_strong = strong_tissue_coloc$locus[strong_tissue_coloc$eqtl_file == "Muscle_Skeletal"]
musk_weak = weak_tissue_coloc$locus[weak_tissue_coloc$eqtl_file == "Muscle_Skeletal"]
loci_musk = c(musk_strong[!musk_strong %in% loci_shared], musk_weak[!musk_weak %in% c(loci_shared, strong_tissue_counts$locus)])

liver_strong = strong_tissue_coloc$locus[strong_tissue_coloc$eqtl_file == "Liver"]
liver_weak = weak_tissue_coloc$locus[weak_tissue_coloc$eqtl_file == "Liver"]
loci_liver = c(liver_strong[!liver_strong %in% loci_shared], liver_weak[!liver_weak %in% c(loci_shared, strong_tissue_counts$locus)])

# Then there's one other special case for adipose...
adp_all = c(adps_strong[adps_strong %in% adpv_strong], adps_weak[(adps_weak %in% adpv_weak) & !(adps_weak %in% strong_tissue_counts$locus)])
loci_adp_both = adp_all[!(adp_all %in% loci_shared)]
# Adjust the other two adps groups to remove this one...
loci_adps = loci_adps[!(loci_adps %in% loci_adp_both)]
loci_adpv = loci_adpv[!(loci_adpv %in% loci_adp_both)]

# Quality check results; make sure everything's been classified
all_loci = c(loci_shared, loci_adps, loci_adpv, loci_musk, loci_liver, loci_adp_both, loci_non_coloc)
stopifnot(dim(coloc_counts)[1] == length(all_loci))
stopifnot(length(table(table(all_loci))) == 1)

step2_list[loci_list %in% loci_adps] = "Adps"
step2_list[loci_list %in% loci_adpv] = "Adpv"
step2_list[loci_list %in% loci_adp_both] = "Adp_Both"
step2_list[loci_list %in% loci_musk] = "Musk"
step2_list[loci_list %in% loci_liver] = "Liver"
step2_list[loci_list %in% loci_shared] = "shared"
step2_list[loci_list %in% loci_non_coloc] = "none"

###################################
# Part 3: Specific to this study
###################################

top_colocs = sub %>% group_by(locus, gwas_trait) %>% summarize(best = max(clpp_mod))
loci_gwas1 = unique(top_colocs[top_colocs$gwas_trait %in% c("ISI_Model_2_AgeSexBMI", "MI_adjBMI") & (top_colocs$best > weak_clpp_threshold),]$locus)

loci_gwas2 = unique(top_colocs[top_colocs$gwas_trait %in% c("FastGlu", "FastInsu_adjBMI", "T2D") & (top_colocs$best > weak_clpp_threshold),]$locus)
loci_gwas2 = loci_gwas2[!loci_gwas2 %in% loci_gwas1]

loci_gwas3 = unique(top_colocs[top_colocs$gwas_trait %in% c("WHRadjBMI") & (top_colocs$best > weak_clpp_threshold),]$locus)
loci_gwas3 = loci_gwas3[!loci_gwas3 %in% c(loci_gwas1, loci_gwas2)]

loci_gwas4 = unique(top_colocs[top_colocs$gwas_trait %in% c("TG", "BMI", "HDL") & (top_colocs$best > weak_clpp_threshold),]$locus)
loci_gwas4 = loci_gwas4[!loci_gwas4 %in% c(loci_gwas1, loci_gwas2, loci_gwas3)]

# These ones are here, but they don't really matter in this analysis because of the whole CHD aspect...will need to
# go back to that original matrix
loci_chd_dont_care = unique(top_colocs[top_colocs$gwas_trait %in% c("CHD") & (top_colocs$best > weak_clpp_threshold),]$locus)
loci_chd_dont_care = loci_chd_dont_care[!loci_chd_dont_care %in% c(loci_gwas1, loci_gwas2, loci_gwas3, loci_gwas4)]

all_loci = c(loci_non_coloc, loci_gwas1, loci_gwas2, loci_gwas3, loci_gwas4, loci_chd_dont_care)
stopifnot(dim(coloc_counts)[1] == length(all_loci))
stopifnot(length(table(table(all_loci))) == 1)

step3_list[loci_list %in% loci_gwas1] = "Tier1"
step3_list[loci_list %in% loci_gwas2] = "Tier2"
step3_list[loci_list %in% loci_gwas3] = "Tier3"
step3_list[loci_list %in% loci_gwas4] = "Tier4"
step3_list[loci_list %in% loci_non_coloc] = "none"

####### Put all the results together ########

classes = data.frame(list(locus=loci_list, step1=step1_list, step2=step2_list, step3=step3_list))
write.table(classes, file="coloc_classification.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

