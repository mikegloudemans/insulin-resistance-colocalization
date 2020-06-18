require(dplyr)

############################################
# Get some stats about individual loci
############################################

# Get full list of results including categories
cats = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", sep="\t", header=TRUE)
snps_per_locus = cats %>% group_by(locus) %>%  summarize(n_snps = length(unique(paste(chr, pos))), named_snps = paste(unique(paste(chr, pos)), collapse=" : "), snp_range = max(pos) - min(pos))

# Plot number of lead snps per locus
# We use half breaks because otherwise 1 and 2 get lumped into the same category
hist(snps_per_locus$n_snps, breaks=seq(0.5,18.5,by=1))
length(snps_per_locus$snp_range[snps_per_locus$n_snps == 1])
length(snps_per_locus$snp_range[snps_per_locus$n_snps == 2])

# Plot histogram of the "widths" of loci
hist(snps_per_locus$snp_range, xlim=c(100000, 5000000), breaks=seq(0,50000000,by=50000))
hist(snps_per_locus$snp_range[snps_per_locus$snp_range != 0], xlim=c(0, 5000000), breaks=seq(0,50000000,by=50000))

# Just for comparison, plot distances between adjacent loci with GWAS SNPs (discarding the ones that cross chromosome boundaries)
locus_windows = cats %>% group_by(chr, locus) %>% summarize(start = min(pos), stop = max(pos))
between_loci = sapply(2:dim(locus_windows)[1], function(x) {locus_windows$start[x] - locus_windows$stop[x-1]})
between_loci = between_loci[between_loci > 0]
hist(between_loci, xlim=c(100000, 5000000), breaks=seq(0,50000000,by=50000))


############################################
# Number of GWAS loci per study / trait
############################################
loci_per_gwas = cats %>% group_by(gwas_short) %>% summarize(num_loci = length(unique(locus)))


############################################
# Get numbers of loci overlapping between
# categories
############################################

# Load info
eqtl_loci = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_only/coloc_classification_2020-05-11.txt", header=TRUE)
sqtl_loci = read.table("output/post_coloc/2020-05-11/refiltered/sqtls_only/coloc_classification_2020-05-11.txt", header=TRUE)
eqtl_or_sqtl_loci = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_classification_2020-05-11.txt", header=TRUE)

# Get number of loci tested for each
dim(eqtl_or_sqtl_loci)
dim(eqtl_loci)
dim(sqtl_loci)

# Number only sQTL loci
dim(eqtl_or_sqtl_loci) - dim(eqtl_loci)
# Only eQTL loci
dim(eqtl_or_sqtl_loci) - dim(sqtl_loci)
# Both (eQTL AND sQTL)
dim(eqtl_loci) - (dim(eqtl_or_sqtl_loci) - dim(sqtl_loci))

# Figure out how many loci fit some specific category
# In this case, which ones have strong colocs of any kind
target_groups = c("loci_0.1", "loci_2", "loci_3")
sum(eqtl_loci$step1 %in% target_groups)
sum(eqtl_loci$step1 %in% target_groups)
sum(eqtl_or_sqtl_loci$step1 %in% target_groups)
# Can then repeat the same math as when we did above to get intersection, etc..

############################################
# How does pancreas factor in to all of this?
############################################
table(eqtl_or_sqtl_loci[c("step2", "pancreas")])
table(eqtl_or_sqtl_loci[c("step1", "pancreas")])
table(eqtl_or_sqtl_loci[c("step3", "pancreas")])

t = table(eqtl_or_sqtl_loci[c("step1", "pancreas")])
strong = t[c(1,4,5),]
rownames(strong) = c("one candidate, one coloc", "multi candidate, one coloc", "multi candidate, multi coloc")
