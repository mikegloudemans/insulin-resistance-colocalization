require(dplyr)

data = read.table("output/post_coloc/de_genes/perturbation_by_coloc.txt", header=TRUE)

# How often is the GWAS directionality not consistent across tissues, for a given gene?
find_inconsistent = data %>% group_by(pert_tissue, perturbation, gene, gwas_file) %>% summarize(num_contexts = length(higher_expression_higher_risk), inconsistent = length(unique(higher_expression_higher_risk))>1)
print(sum(find_inconsistent$inconsistent))
# At least in current form, it's never inconsistent

# What about even across GWAS? Do we ever see inconsistency across different GWAS?
find_inconsistent_aggregated = data %>% group_by(pert_tissue, perturbation, gene) %>% summarize(num_contexts = length(higher_expression_higher_risk), inconsistent = length(unique(higher_expression_higher_risk))>1)
print(sum(find_inconsistent_aggregated$inconsistent))
inconsistent_aggregated_genes = find_inconsistent_aggregated %>% filter(inconsistent==TRUE)
# When there's inconsistency but it's a different lead variant, this could be an LD tagging issue
# rather than actual directional discordance, to be fair

# TODO: Are there any GWAS that are seeming inconsistent more often than others?
# TODO: How often are each group of pairwise GWAS directionally concordant when appearing with same gene?

# TODO: For each tissue / perturbation / GWAS combo, plot the percentage of all DE genes that were
# 	DE in the risk-increasing direction (heatmap works for this, can be colored by direction)
per_gene_directions = data %>% group_by(pert_tissue, perturbation, gene, gwas_file) %>% summarize(pert_increases_risk=sum(perturbation_increases_risk=="True") > sum(perturbation_increases_risk=="False"), pert_decreases_risk = sum(perturbation_increases_risk=="False") > sum(perturbation_increases_risk=="True"), expression_increases_risk=sum(higher_expression_higher_risk=="True") > sum(higher_expression_higher_risk=="False"), expression_decreases_risk = sum(higher_expression_higher_risk=="False") > sum(higher_expression_higher_risk=="True"))

# Do this in a more thorough way, since it's potentially interesting...
sum(per_gene_directions$pert_increases_risk)
sum(per_gene_directions$pert_decreases_risk)
sum(per_gene_directions$expression_increases_risk)
sum(per_gene_directions$expression_decreases_risk)

expression_effects = data %>% group_by(gwas_file, gene) %>% summarize(expression_increases_risk=sum(higher_expression_higher_risk=="True") > sum(higher_expression_higher_risk=="False"), expression_decreases_risk = sum(higher_expression_higher_risk=="False") > sum(higher_expression_higher_risk=="True"))
per_gwas_expression_effects = expression_effects %>% group_by(gwas_file) %>% summarize(pos_pos = sum(expression_increases_risk), pos_neg = sum(expression_decreases_risk))

#frac_positive_per_perturbation = per_gene_directions %>% group_by(pert_tissue, perturbation, gwas_file) %>% summarize(fraction = sum(pos_direction) / sum(pos_direction + neg_direction))

# How often is the tissue the same for pert and coloc? I don't know...

###########################################################################

# DONE: Double-check GWAS effect directions

# - Thought: should we re-call DE genes from the original set using only the ones we're interested in here?
#            I think this would be legit, though not sure how much it would increase the number we're including then
# - Should we include loci with multiple colocs? Maybe take just the top gene from each though?
# 	- Or only loci with a certain number of colocs
# - Choose the relevant tissue / trait pairs, and only focus on those...including filter to the relevant tissue colocalizations

# Tile plots for single-coloc genes
#	- Supplement. + correct for lower FDR (supplement w/ all overlaps)
#	- Could subset down to only a few
#		- Important genes
#		- Genes relevant with other perturbations

##############################################

# Probably not for this:
# - Could pick top CLPP variant instead of top GWAS variant

# TODO? Probably not
# For each tissue / perturbation combo
# For each perturbation / GWAS combo (across all tissues)
# For each perturbation (across all tissues)
