require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

#####################################################
### Setup and loading config settings
#####################################################

# Min distance between two adjacent GWAS SNPs to be considered different loci
# Can be overridden in config file
default_min_locus_distance = 1000000

config_file = commandArgs(trailingOnly=TRUE)[1]

# Load pre-specified config file
config = fromJSON(file=config_file)

if ("min_locus_distance" %in% names(config))
{
	min_locus_distance = as.numeric(config$min_locus_distance)
} else
{
	min_locus_distance = default_min_locus_distance
}

# In this file, we'll store some numbers summarizing the number
# of loci remaining after each filtering step
# (This may eventually be part of a figure)
write_summary = function(x, append=TRUE)
{
	summary_file = paste0(config$out_dir, "/", "filtering_summary_counts.txt")
	write(x, file=summary_file, append=append)
}
write_summary("Summary of filters:\n--------------------\n", append=FALSE)

# We select the traits we are interested in. We originally load
# all traits but then subset down to these.
trait_list = sapply(config$kept_traits, function(x) {s=strsplit(x, "/"); return(s[[1]][length(s[[1]])])})
names(trait_list) = NULL

# create a function that tests GWAS pval / eQTL trait rows for validity
pval_passing = function(x, threshold_set)
{
	pvalue = as.numeric(x[1])
	trait = unlist(x[2])
	if (pvalue < threshold_set$standard)
	{
		return(TRUE)
	}
	else if (("exceptions" %in% names(threshold_set)) && (trait %in% names(threshold_set$exceptions)))
	{
		if (pvalue < threshold_set$exceptions[[trait]])
		{
			return(TRUE)
		}
	}
	return(FALSE)
}

print("Traits included")
print(trait_list)

#####################################################
### Summarize pre-colocalization filtering
#####################################################

# Count the total number of unique SNPs being tested that were the lead SNP for at least one GWAS.
# Note: We filter out the SNPs that were selected on the basis of excluded GWAS, meaning the only
# ones we have left will be significant in at least one of the GWAS we care about (that was an
# inclusion criterion for this list).

### Get number of unique GWAS lead SNPs.
all_gwas_snps = read.table(config$gwas_snp_list, col.names=c("chr", "snp.pos", "pvalue", "source_trait", "source_file"), fill=TRUE, skip=1)

## Temporary, for local runs ##
all_gwas_snps = all_gwas_snps[all_gwas_snps$source_file %in% config$kept_traits,]
all_gwas_snps = all_gwas_snps[apply(all_gwas_snps[c("pvalue", "source_trait")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]
all_gwas_snps$ref_snp = paste(all_gwas_snps$chr, all_gwas_snps$snp.pos, sep="_")

write_summary(sprintf("Unique GWAS lead SNPs:\t%d\n", length(unique(all_gwas_snps$ref_snp))))


### Count the number of SNPs considered for each GWAS.

hits_per_trait = all_gwas_snps %>% group_by(source_trait) %>% summarize(snp_count = length(unique(ref_snp)))
hits_per_trait$trait = sapply(as.character(hits_per_trait$source_trait), function(x) {s=strsplit(x, "/"); return(s[[1]][length(s[[1]])])})
print(hits_per_trait)

# The new algorithm for grouping to loci!
# Invariant to the total number of variants in the list

group_to_loci = function(x)
{
	# Get the set of unique reference SNPs
	ids = unique(as.character(x))
	ids = ids[order(ids)]
        chr = as.numeric(sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][1]}))
        pos = as.numeric(sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][2]}))
        loc_nums = rep(0, length(ids))

	# I ran liftOver to convert to hg38
	# liftOver fourier_ls-all.hg19.bed /mnt/lab_data/montgomery/shared/liftOver/chains/hg19ToHg38.over.chain.gz fourier_ls-all.hg38.bed fourier_ls-all.hg38.failed.bed 

	# Load European independent LD block partitioning from LDetect
	ldetect = read.table("data/ldetect/fourier_ls-all.hg38.connected.bed", header=FALSE)
	colnames(ldetect) = c("chr", "start", "stop")
	ldetect$chr = as.numeric(gsub("chr", "", ldetect$chr))
	ldetect$locus = 1:dim(ldetect)[1]

	# Assign each SNP to its own locus segment
        for (i in 1:length(ids))
	{
		# Get the chromosome we're looking for
		sub_chr = ldetect[ldetect$chr == chr[i],]
		# Get the exact block we're looking for
		sub_locus = sub_chr[sub_chr$start <= pos[i] & sub_chr$stop > pos[i],]
		loc_nums[i] = sub_locus$locus[1]
	}

	# TODO: Deal with any NA's
	mapped_loci = sapply(x, function(j) 
	{
			loc_nums[which(ids == j)]
        })

	return(mapped_loci)
}

### How many individual loci are represented in our collection of SNP?

all_gwas_snps$locus = group_to_loci(all_gwas_snps$ref_snp)


##############################################################################
### Count the total number of genes found before filtering by eQTL p-values.
##############################################################################

# NOTE: Again, to be in this list a SNP-gene pair must be significant at least in the GWAS,
# so no SNPs that shouldn't be tested will slip through this step.

all_snp_gene_pairs = read.table(config$snp_gene_pair_list, col.names=c("chr", "snp.pos", "source_pvalue", "source_trait", "lookup_pvalue", "source_file", "feature", "lookup_file"), skip=1)

# A bit of stuff just for more understandable labels
all_snp_gene_pairs$ensembl = all_snp_gene_pairs$feature
all_snp_gene_pairs$ensembl = gsub(":", ".", all_snp_gene_pairs$ensembl)

# Map splice events to the corresponding gene
splice_map = read.table("data/sqtls/gtex_v8/splice_to_gene_map.txt.gz", header=TRUE)
splice_map$splice_junction = gsub(":", ".", splice_map$splice_junction)
splice_indices = grepl("clu", all_snp_gene_pairs$feature)
all_snp_gene_pairs$ensembl[splice_indices] = as.character(splice_map$gene[match(all_snp_gene_pairs$ensembl[splice_indices], splice_map$splice_junction)])

# Trim the ensemble version number
all_snp_gene_pairs$ensembl = sapply(as.character(all_snp_gene_pairs$ensembl), function(x) {strsplit(x, "\\.")[[1]][1]})

## Temporary, for local runs ##
#all_snp_gene_pairs$gwas_file = gsub("/users/mgloud/projects/insulin_resistance/", "", all_snp_gene_pairs$gwas_file)
#all_snp_gene_pairs$trait = gsub("/users/mgloud/projects/insulin_resistance/", "", all_snp_gene_pairs$trait)
#all_snp_gene_pairs$eqtl_file = gsub("/users/mgloud/projects/brain_gwas/", "", all_snp_gene_pairs$eqtl_file)

all_snp_gene_pairs = all_snp_gene_pairs[all_snp_gene_pairs$source_file %in% config$kept_traits,]
all_snp_gene_pairs = all_snp_gene_pairs[all_snp_gene_pairs$lookup_file %in% config$kept_eqtls,]
all_snp_gene_pairs$ref_snp = paste(all_snp_gene_pairs$chr, all_snp_gene_pairs$snp.pos, sep="_")
all_snp_gene_pairs = all_snp_gene_pairs[apply(all_snp_gene_pairs[c("source_pvalue", "source_trait")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]
all_snp_gene_pairs = all_snp_gene_pairs[all_snp_gene_pairs$ref_snp %in% unique(all_gwas_snps$ref_snp),]
all_snp_gene_pairs$locus = group_to_loci(all_snp_gene_pairs$ref_snp)


##############################################################################
### Count the list of tests we actually tried to run (after eQTL thresholding).
##############################################################################

# This also tells us which ones we've lost.
# Once more, after filtering, SNPs will only be in this list if they were significant in
# at least one of the GWAS we care about.

testable_snps = read.table(config$coloc_test_list, header=TRUE)

testable_snps$source_file = testable_snps$source_file
testable_snps$lookup_file = testable_snps$lookup_file
testable_snps$source_file = gsub("/users/mgloud/projects/insulin_resistance/", "", testable_snps$source_file)
testable_snps$source_trait = gsub("/users/mgloud/projects/insulin_resistance/", "", testable_snps$source_trait)
testable_snps$lookup_file = gsub("/users/mgloud/projects/brain_gwas/", "", testable_snps$lookup_file)

testable_snps = testable_snps[testable_snps$source_file %in% config$kept_traits,]
testable_snps = testable_snps[testable_snps$lookup_file %in% config$kept_eqtls,]
print(dim(testable_snps))

# A bit of stuff just for more understandable labels
testable_snps$ensembl = testable_snps$lookup_trait
testable_snps$ensembl = gsub(":", ".", testable_snps$ensembl)

# Map splice events to the corresponding gene
splice_indices = grepl("clu", testable_snps$lookup_trait)
testable_snps$ensembl[splice_indices] = as.character(splice_map$gene[match(testable_snps$ensembl[splice_indices], splice_map$splice_junction)])

# Trim the ensemble version number
testable_snps$ensembl = sapply(as.character(testable_snps$ensembl), function(x) {strsplit(x, "\\.")[[1]][1]})

testable_snps$ref_snp = paste(testable_snps$chr, testable_snps$snp_pos, sep="_")
testable_snps = testable_snps[apply(testable_snps[c("source_pvalue", "source_trait")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]
testable_snps = testable_snps[apply(testable_snps[c("lookup_pvalue", "lookup_file")], 1, FUN=pval_passing, threshold = config$eqtl_pval_threshold),]
testable_snps = testable_snps[testable_snps$ref_snp %in% unique(all_gwas_snps$ref_snp),]
testable_snps$locus = group_to_loci(testable_snps$ref_snp)
testable_snps$source_trait = sapply(as.character(testable_snps$source_trait), function(x) {s=strsplit(x, "/"); return(s[[1]][length(s[[1]])])})
print(dim(testable_snps))

### How many SNPs per trait are now remaining?

remaining_hits_per_trait = testable_snps %>% group_by(source_trait) %>% summarize(snp_count = length(unique(ref_snp)))
print(data.frame(remaining_hits_per_trait))

# Figure out which SNPs we dropped from the first list by eQTL thresholding
dropped_snps = all_gwas_snps[!(all_gwas_snps$ref_snp %in% unique(testable_snps$ref_snp)),]

# How many SNPs dropped per trait after applying eQTL filters?
dropped_hits_per_trait = dropped_snps %>% group_by(source_trait) %>% summarize(snp_count = length(unique(ref_snp)))
print(dropped_hits_per_trait)

############################################################
### Post-colocalization summarization and QC
############################################################

# We load all the colocalization files.

tabs = list()
err_tab = list()
skip_tab = list()
i = 1

folders = config$coloc_out_dirs

for (folder in folders)
{
	files = dir(folder)
	error = files[grep("ERROR", files)]
	skip = files[grep("skip", files)]
	files = files[grep("clpp" ,files)]
	
	if (length(error) > 0)
	{
		err_tab[[i]] = read.table(paste(folder, error, sep="/"), header=FALSE, sep="\t", fill=TRUE, col.names = c("gwas_file", "eqtl_file", "snp.chrom", "snp.pos", "restrict_gene", "trait", "error", "error2"))
		err_tab[[i]]$error = paste(err_tab[[i]]$error, err_tab[[i]]$error2, sep="|")
	}
	if (length(skip) > 0)
	{
		skip_tab[[i]] = read.table(paste(folder, skip, sep="/"), header=FALSE, sep="\t", fill=TRUE, col.names = c("gwas_file", "eqtl_file", "snp.chrom", "snp.pos", "feature", "error", "gwas_data"))
	}

	for (file in files)
	{
		tabs[[i]] = read.table(paste(folder, file, sep="/"), header=FALSE, col.names=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "n_snps", "clpp", "neg_log_gwas_pval","neg_log_eqtl_pval","base_gwas_file","clpp_mod"))
		i = i + 1
	}
}

results = do.call(rbind, tabs)
errors = do.call(rbind, err_tab)
skips = do.call(rbind, skip_tab)

# It's possible will have exact duplicates if the traits weren't run all at once.
# I fix this here.
results = results[!duplicated(results),]

# Rename GWAS and eQTL files for complete consistency
# I don't like this solution for the long term, but
# since we're dealing with files that don't have the info we
# want, this is probably the best we can do right now.
results$base_gwas_file = as.character(results$base_gwas_file)
results$eqtl_file = as.character(results$eqtl_file)

#results$base_gwas_file = gsub("/users/mgloud/projects/insulin_resistance/", "", results$base_gwas_file)
#results$eqtl_file = gsub("/users/mgloud/projects/brain_gwas/", "", results$eqtl_file)
#errors$gwas_file = gsub("/users/mgloud/projects/insulin_resistance/", "", errors$gwas_file)
#errors$eqtl_file = gsub("/users/mgloud/projects/brain_gwas/", "", errors$eqtl_file)
#skips$gwas_file = gsub("/users/mgloud/projects/insulin_resistance/", "", skips$gwas_file)
#skips$eqtl_file = gsub("/users/mgloud/projects/brain_gwas/", "", skips$eqtl_file)

for (alias in names(config$gwas_aliases))
{
	if (alias %in% (unique(results$base_gwas_file)))
	{
		results[results$base_gwas_file == alias,]$base_gwas_file = config$gwas_aliases[[alias]]
	}
}
for (alias in names(config$eqtl_aliases))
{
	if (alias %in% (unique(results$eqtl_file)))
	{
		results[results$eqtl_file == alias,]$eqtl_file = config$eqtl_aliases[[alias]]
	}
}

# Filter down to the ones that are in kept traits list
results = results[results$base_gwas_file %in% config$kept_traits,]
results = results[results$eqtl_file %in% config$kept_eqtls,]
errors = errors[errors$gwas_file %in% config$kept_traits,]
errors = errors[errors$eqtl_file %in% config$kept_eqtls,]
skips = skips[skips$gwas_file %in% config$kept_traits,]
skips = skips[skips$eqtl_file %in% config$kept_eqtls,]

skips$ref_snp = paste(skips$snp.chrom, skips$snp.pos, sep="_")
errors$ref_snp = paste(errors$snp.chrom, errors$snp.pos, sep="_")

# A bit of stuff just for more understandable labels
results$ensembl = as.character(results$feature)

# Map splice events to the corresponding gene
splice_indices = grepl("clu", results$ensembl)
results$ensembl[splice_indices] = as.character(splice_map$gene[match(results$ensembl[splice_indices], splice_map$splice_junction)])

# Trim the ensemble version number
results$ensembl = sapply(as.character(results$ensembl), function(x) {strsplit(x, "\\.")[[1]][1]})

results$gwas_short = "none"
results$eqtl_short = "none"
if ("gwas_short_names" %in% names(config))
{
	for (full_name in names(config$gwas_short_names))
	{
		results$gwas_short[results$base_gwas_file == full_name] = config$gwas_short_names[[full_name]]
	}
}
if ("eqtl_short_names" %in% names(config))
{
	for (full_name in names(config$eqtl_short_names))
	{
		results$eqtl_short[results$eqtl_file == full_name] = config$eqtl_short_names[[full_name]]
	}
}

results$gwas_pval = 10^(-as.numeric(results$neg_log_gwas_pval))
results$eqtl_pval = 10^(-as.numeric(results$neg_log_eqtl_pval))

# These ones are dropped because they may have been included at first, but
# didn't pass our final GWAS / eQTL significance thresholds.
results = results[apply(results[c("gwas_pval", "base_gwas_file")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),]
results = results[apply(results[c("eqtl_pval", "eqtl_file")], 1, FUN=pval_passing, threshold = config$eqtl_pval_threshold),]
results = results[results$ref_snp %in% unique(all_gwas_snps$ref_snp),]
print("dim(results)")
dim(results)

results$locus = group_to_loci(results$ref_snp)

coloc_passing = results[as.numeric(results$clpp_mod) >= as.numeric(config$target_clpp_mod_cutoff),]
coloc_passing = results[as.numeric(results$n_snps) >= as.numeric(config$target_num_snps),]

#########################################################
## We need to filter all lists to make sure they've
## removed the loci that were dropped for coloc analysis
######################################################### 

# This solves the problem
testable_snps = testable_snps[testable_snps$ref_snp %in% unique(results$ref_snp),]

#########################################################
## Now we write summaries to a file
#########################################################

write_summary(sprintf("Unique GWAS loci:\t%d\n", length(unique(all_gwas_snps$locus))))
write_summary(sprintf("Total number of unique pairs of SNP and GWAS trait\t%d\n", length(unique(paste(all_gwas_snps$ref_snp, all_gwas_snps$source_trait, sep="_")))))
write_summary(sprintf("Total number of unique pairs of locus and GWAS trait\t%d\n", length(unique(paste(all_gwas_snps$locus, all_gwas_snps$source_trait, sep="_")))))

write_summary(sprintf("Total number of unique features overlapping GWAS\t%d\n", length(unique(all_snp_gene_pairs$feature))))
write_summary(sprintf("Total number of unique genes overlapping GWAS\t%d\n", length(unique(all_snp_gene_pairs$ensembl))))
write_summary(sprintf("Total number of unique SNP-gene pairs\t%d\n", length(unique(paste(all_snp_gene_pairs$ref_snp, all_snp_gene_pairs$ensembl, sep="_")))))
write_summary(sprintf("Total number of unique locus-gene pairs\t%d\n", length(unique(paste(all_snp_gene_pairs$locus, all_snp_gene_pairs$ensembl, sep="_")))))

write_summary(sprintf("Number of testable SNPs after filtering with eQTLs:\t%d\n", length(unique(testable_snps$ref_snp))))
write_summary(sprintf("Number of testable loci after filtering:\t%d\n", length(unique(testable_snps$locus))))
write_summary(sprintf("Number of testable features after filtering:\t%d\n", length(unique(testable_snps$lookup_trait))))
write_summary(sprintf("Number of testable genes after filtering:\t%d\n", length(unique(testable_snps$ensembl))))
write_summary(sprintf("Number of unique SNP-gene pairs to test after filtering:\t%d\n", length(unique(paste(testable_snps$ref_snp, testable_snps$ensembl, sep="_")))))
write_summary(sprintf("Number of unique locus-gene pairs to test after filtering:\t%d\n", length(unique(paste(testable_snps$locus, testable_snps$ensembl, sep="_")))))

write_summary(sprintf("Number of colocalized lead SNPs:\t%d\n", length(unique(coloc_passing$ref_snp))))
write_summary(sprintf("Number of colocalized loci:\t%d\n", length(unique(coloc_passing$locus))))
write_summary(sprintf("Number of colocalized features:\t%d\n", length(unique(coloc_passing$feature))))
write_summary(sprintf("Number of colocalized genes:\t%d\n", length(unique(coloc_passing$ensembl))))
write_summary(sprintf("Number of colocalized SNP-gene pairs:\t%d\n", length(unique(paste(coloc_passing$ref_snp, coloc_passing$ensembl, sep="_")))))
write_summary(sprintf("Number of colocalized locus-gene pairs:\t%d\n", length(unique(paste(coloc_passing$locus, coloc_passing$ensembl, sep="_")))))

write.table(all_gwas_snps, file=paste(config$out_dir, "all_gwas_snps.txt", sep="/"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
write.table(all_snp_gene_pairs, file=paste(config$out_dir, "all_snp_gene_pairs_prefiltering.txt", sep="/"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# TODO: Note this part is really slow, might want to come up with a faster way of achieving this
# I think though, it's just because there are a ton of files to look through, so not much way around that.

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

write.table(results, file=paste(config$out_dir, "full_coloc_results_qced.txt", sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
