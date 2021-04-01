suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(require(rjson)))

coloc_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt"

data = read.table(coloc_file, header=TRUE, sep="\t", fill=TRUE)
eqtls = data[grepl("eQTLs", data$eqtl_file),]
sqtls = data[grepl("sQTLs", data$eqtl_file),]

###################################################
# Part 1: eCDFs of CLPP scores and CLPP mod
###################################################

# For both CLPP and CLPP-mod...

# CLPP scores CDF
plot(ecdf(data$clpp_mod), xlab="CLPP mod", ylab="quantile", main="CLPP mod eCDF", lwd=2)
for(i in seq(0,1,0.05))
{
	abline(lty="dotted", col="lightgray", v=i)
	abline(lty="dotted", col="lightgray", a=i, b=0)
}
for(i in seq(0,1,0.2))
{
	abline(col="lightgray", v=i)
	abline(col="lightgray", a=i, b=0)
}

# CLPP scores per GWAS CDFs?
plot(ecdf(data$clpp), xlab="CLPP", ylab="quantile", main="CLPP (normal) eCDF (truncated x-axis)", lwd=2, xlim=c(0,0.05))
plot(ecdf(data$clpp), xlab="CLPP", ylab="quantile", main="CLPP (normal) eCDF (super truncated x-axis)", lwd=2, xlim=c(0,0.01))
plot(ecdf(data$clpp), xlab="CLPP", ylab="quantile", main="CLPP (normal) eCDF (super super truncated x-axis)", lwd=2, xlim=c(0,0.001))

# CLPP scores per GWAS CDFs?
plot(ecdf(data$clpp_mod), xlab="CLPP mod", ylab="quantile", main="CLPP mod eCDF, split by GWAS", lwd=2)
for (gwas in unique(data$gwas_short))
{
	sub = data[data$gwas_short == gwas,]
	lines(ecdf(sub$clpp_mod), lwd=1, cex.points=0)
}

# CLPP scores per eQTL CDFs?
plot(ecdf(data$clpp_mod), xlab="CLPP mod", ylab="quantile", main="CLPP mod eCDF, split by QTL tissue and type", lwd=2)
lines(ecdf(sqtls$clpp_mod), lwd=2, col="red")
lines(ecdf(eqtls$clpp_mod), lwd=2, col="blue")
for (eqtl in unique(sqtls$eqtl_short))
{
	sub = sqtls[sqtls$eqtl_short == eqtl,]
	lines(ecdf(sub$clpp_mod), lwd=1, cex.points=0, col="red")
}
for (eqtl in unique(eqtls$eqtl_short))
{
	sub = eqtls[eqtls$eqtl_short == eqtl,]
	lines(ecdf(sub$clpp_mod), lwd=1, cex.points=0, col="blue")
}

###################################################
# Part 2: eCDFs of loci and genes detected at
# various cutoffs of CLPP scores and CLPP mod
###################################################

cutoffs = seq(0, 1, by=0.001)

# Gene CDFs

# CLPP mod

genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp_mod >= x,]$ensembl))})
sqtl_genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(sqtls[sqtls$clpp_mod >= x,]$ensembl))})
eqtl_genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(eqtls[eqtls$clpp_mod >= x,]$ensembl))})

plot(cutoffs, genes_per_cutoff, type="l", xlab="CLPP mod cutoff", ylab="# unique colocalized genes", main="Varying CLPP mod thresholds")

for(i in seq(0,1,0.05))
{
	abline(lty="dotted", col="lightgray", v=i)
}
for(i in seq(0,1,0.2))
{
	abline(col="lightgray", v=i)
}
for(i in seq(0, 4000, 100))
{
	abline(lty="dotted", col="lightgray", a=i, b=0)
}
for(i in seq(0, 4000 , 500))
{
	abline(col="lightgray", a=i, b=0)
}
lines(cutoffs, genes_per_cutoff)
lines(cutoffs, sqtl_genes_per_cutoff, col="red")
lines(cutoffs, eqtl_genes_per_cutoff, col="blue")

# CLPP

genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp >= x,]$ensembl))})
sqtl_genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(sqtls[sqtls$clpp >= x,]$ensembl))})
eqtl_genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(eqtls[eqtls$clpp >= x,]$ensembl))})

plot(cutoffs, genes_per_cutoff, type="l", xlab="CLPP cutoff", ylab="# unique colocalized genes", main="Varying CLPP thresholds", xlim=c(0,0.2), ylim=c(0,1500))
for(i in seq(0,0.2,0.01))
{
	abline(lty="dotted", col="lightgray", v=i)
}
for(i in seq(0,0.2,0.05))
{
	abline(col="lightgray", v=i)
}
for(i in seq(0, 2000, 100))
{
	abline(lty="dotted", col="lightgray", a=i, b=0)
}
for(i in seq(0, 2000 , 500))
{
	abline(col="lightgray", a=i, b=0)
}

lines(cutoffs, genes_per_cutoff)
lines(cutoffs, sqtl_genes_per_cutoff, col="red")
lines(cutoffs, eqtl_genes_per_cutoff, col="blue")


# Locus CDFs

loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp_mod >= x,]$locus))})
sqtl_loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(sqtls[sqtls$clpp_mod >= x,]$locus))})
eqtl_loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(eqtls[eqtls$clpp_mod >= x,]$locus))})

plot(cutoffs, loci_per_cutoff, type="l", xlab="CLPP mod cutoff", ylab="# unique colocalized loci", main="Varying CLPP mod thresholds")

for(i in seq(0,1,0.05))
{
	abline(lty="dotted", col="lightgray", v=i)
}
for(i in seq(0,1,0.2))
{
	abline(col="lightgray", v=i)
}
for(i in seq(0, 800, 10))
{
	abline(lty="dotted", col="lightgray", a=i, b=0)
}
for(i in seq(0, 800, 100))
{
	abline(col="lightgray", a=i, b=0)
}
lines(cutoffs, loci_per_cutoff)
lines(cutoffs, sqtl_loci_per_cutoff, col="red")
lines(cutoffs, eqtl_loci_per_cutoff, col="blue")

# CLPP

loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp >= x,]$locus))})
sqtl_loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(sqtls[sqtls$clpp >= x,]$locus))})
eqtl_loci_per_cutoff = sapply(cutoffs, function(x) {length(unique(eqtls[eqtls$clpp >= x,]$locus))})

plot(cutoffs, loci_per_cutoff, type="l", xlab="CLPP cutoff", ylab="# unique colocalized loci", main="Varying CLPP thresholds", xlim=c(0,0.2))
for(i in seq(0,0.2,0.01))
{
	abline(lty="dotted", col="lightgray", v=i)
}
for(i in seq(0,0.2,0.05))
{
	abline(col="lightgray", v=i)
}
for(i in seq(0, 100, 10))
{
	abline(lty="dotted", col="lightgray", a=i, b=0)
}
for(i in seq(0, 800 , 100))
{
	abline(col="lightgray", a=i, b=0)
}

lines(cutoffs, loci_per_cutoff)
lines(cutoffs, sqtl_loci_per_cutoff, col="red")
lines(cutoffs, eqtl_loci_per_cutoff, col="blue")

########################################################
# Part 3: Direct comparison of CLPP and CLPP mod
########################################################

# CLPP mod vs. CLPP scatterplot
plot(data$clpp_mod, data$clpp, pch=18)
cor.test(data$clpp_mod, data$clpp)
plot(rank(data$clpp_mod), rank(data$clpp), pch=18)
cor.test(rank(data$clpp_mod), rank(data$clpp))

# CLPP mod and CLPP vs. number SNPs tested at locus
# The plots are a mess though
#plot(rank(data$clpp), rank(data$n_snps), pch=18)
cor.test(rank(data$clpp), rank(data$n_snps))
#plot(rank(data$clpp_mod), rank(data$n_snps), pch=18)
cor.test(rank(data$clpp_mod), rank(data$n_snps))

# vs. GWAS pval
plot(rank(data$clpp), rank(data$min_gwas_pval), pch=18)
cor.test(rank(data$clpp), rank(data$min_gwas_pval))
plot(rank(data$clpp_mod), rank(data$min_gwas_pval), pch=18)
cor.test(rank(data$clpp_mod), rank(data$min_gwas_pval))

# vs. eQTL pval
plot(rank(data$clpp), rank(data$min_eqtl_pval), pch=18)
cor.test(rank(data$clpp), rank(data$min_eqtl_pval))
plot(rank(data$clpp_mod), rank(data$min_eqtl_pval), pch=18)
cor.test(rank(data$clpp_mod), rank(data$min_eqtl_pval))

# vs max (GWAS, eQTL pval)
plot(rank(data$clpp), rank(pmin(data$min_eqtl_pval, data$min_gwas_pval)), pch=18)
cor.test(rank(data$clpp), rank(pmin(data$min_eqtl_pval, data$min_gwas_pval)))
plot(rank(data$clpp_mod), rank(pmin(data$min_eqtl_pval, data$min_gwas_pval)), pch=18)
cor.test(rank(data$clpp_mod), rank(pmin(data$min_eqtl_pval, data$min_gwas_pval)))

########################################################
# Part 4: Direct comparison of CLPP and CLPP mod
# specifically for LD patterns
########################################################

# NOTE: This might not be 100% reliable when the lead variant isn't actually
# in the LD lookup table to begin with -- will work for some at least though
#
# Also it matters more how many of these variants were actually tested, not
# the absolute total number of them
#

# CLPP mod and CLPP vs. number LD buddies
ld_file = "data/ld/EUR_geno_hg38.txt.gz"

data$ld_buddy_count = -1

count_ld_buddies = function(chr, snp_pos)
{
		# Get all SNPs that are LD buddies of any of those lead SNPs (r2>0.8) and add them to a "set" for the locus
		ld = length(unlist(system(sprintf("tabix %s chr%s:%s-%s", ld_file, as.character(chr), 
					   as.character(snp_pos), as.character(snp_pos)), intern=TRUE)))
		return(ld)
}

all_snps = unique(data[c("chr", "pos")])
for (i in 1:dim(all_snps)[1])
{
	chr = all_snps$chr[i]
	pos = all_snps$pos[i]
	count = count_ld_buddies(chr, pos)

	data$ld_buddy_count[(data$chr = chr) & (data$pos == pos)] = count
}

plot(data$ld_buddy_count, data$clpp, pch=18)
cor.test(data$ld_buddy_count, data$clpp)
plot(rank(data$ld_buddy_count), rank(data$clpp), pch=18)
cor.test(rank(data$ld_buddy_count), rank(data$clpp))

plot(data$ld_buddy_count, data$clpp_mod, pch=18)
cor.test(data$ld_buddy_count, data$clpp_mod)
plot(rank(data$ld_buddy_count), rank(data$clpp_mod), pch=18)
cor.test(rank(data$ld_buddy_count), rank(data$clpp_mod))

if (FALSE)
{
	# Worth checking, but doesn't really make a big difference

	sub = data[data$ld_buddy_count != 0,]

	plot(sub$ld_buddy_count, sub$clpp, pch=18)
	cor.test(sub$ld_buddy_count, sub$clpp)
	plot(rank(sub$ld_buddy_count), rank(sub$clpp), pch=18)
	cor.test(rank(sub$ld_buddy_count), rank(sub$clpp))

	plot(sub$ld_buddy_count, sub$clpp_mod, pch=18)
	cor.test(sub$ld_buddy_count, sub$clpp)
	plot(rank(sub$ld_buddy_count), rank(sub$clpp), pch=18)
	cor.test(rank(sub$ld_buddy_count), rank(sub$clpp))
}

# Let's maybe just get one score per locus, for simplicity
locus_groups = data %>% group_by(chr, pos) %>% summarize(ld_buddy_count = median(ld_buddy_count), best_clpp = max(clpp), best_clpp_mod = max(clpp_mod))

hist(locus_groups$ld_buddy_count, breaks=seq(0,max(locus_groups$ld_buddy_count)+1, 1)-0.1, xlim=c(0, 100))

plot(rank(locus_groups$ld_buddy_count), rank(locus_groups$best_clpp), pch=18)
cor.test(rank(locus_groups$ld_buddy_count), rank(locus_groups$best_clpp))
# It kind of seems like even if there's not a great correlation overall, the ones with some colocalization tend to be sensitive to these numbers
# So it might actually have an effect depending on where the threshold is set...
some_locus_groups = locus_groups[locus_groups$best_clpp > 0.01,]
plot(rank(some_locus_groups$ld_buddy_count), rank(some_locus_groups$best_clpp), pch=18)
cor.test(rank(some_locus_groups$ld_buddy_count), rank(some_locus_groups$best_clpp))

plot(rank(locus_groups$ld_buddy_count), rank(locus_groups$best_clpp_mod), pch=18)
cor.test(rank(locus_groups$ld_buddy_count), rank(locus_groups$best_clpp_mod))
# It kind of seems like even if there's not a great correlation overall, the ones with some colocalization tend to be sensitive to these numbers
# So it might actually have an effect depending on where the threshold is set...
some_locus_groups = locus_groups[locus_groups$best_clpp_mod > 0.3,]
plot(rank(some_locus_groups$ld_buddy_count), rank(some_locus_groups$best_clpp_mod), pch=18)
cor.test(rank(some_locus_groups$ld_buddy_count), rank(some_locus_groups$best_clpp_mod))


########################################################
# Part 5: Messing around with FDRs
########################################################

# We can add this to the CDFs actually, it's just cumsum(1-prob I think)
# And can additionally then make a plot showing it as percents...

cutoffs = seq(0, 1, by=0.001)

# NOTE: Hard to estimate the exact number of GENES that will be falsely discovered,
# we can really only estimate the exact number of coloc events that will be false.

# With CLPP mod
# How about the number of coloc events total though?
number_results = sapply(cutoffs, function(x) {sum(data$clpp_mod >= x)})
expected_false_positives = ceiling(sapply(cutoffs, function(x) {sum(1-(data[data$clpp_mod >= x,]$clpp_mod))}))
expected_fpr = expected_false_positives / number_results

genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp_mod >= x,]$ensembl))})
plot(cutoffs, genes_per_cutoff, type="l", xlab="CLPP mod cutoff", ylab="# unique colocalized genes", main="Varying CLPP mod thresholds")
lines(cutoffs, expected_fpr * 4000, lty="dashed")


# With old CLPP
number_results = sapply(cutoffs, function(x) {sum(data$clpp >= x)})
expected_false_positives = ceiling(sapply(cutoffs, function(x) {sum(1-(data[data$clpp >= x,]$clpp))}))
expected_fpr = expected_false_positives / number_results

genes_per_cutoff = sapply(cutoffs, function(x) {length(unique(data[data$clpp >= x,]$ensembl))})
plot(cutoffs, genes_per_cutoff, type="l", xlab="CLPP cutoff", ylab="# unique colocalized genes", main="Varying CLPP thresholds", xlim=c(0,0.2), ylim=c(0,1500))
lines(cutoffs, expected_fpr * 1500, lty="dashed")


