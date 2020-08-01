require(dplyr)

# Take a list of known monogenic genes and see if they were
# tested for colocalization, and if they colocalized.

# Could also see if they were even in GTEx annotations,
# if they were expressed, where they fell out, etc,
# but not going to worry about that right this moment.

monogenic = read.table("data/monogenic_ir_gene_list.txt", header=FALSE, sep="\t")
colnames(monogenic) = c("HGNC", "ensembl")
monogenic = monogenic[!duplicated(monogenic),]

coloc = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", sep="\t", header=TRUE)

# How many of the monogenic genes were or weren't tested for colocalization
# in our analysis?
sum(monogenic$ensembl %in% coloc$ensembl)
sum(!(monogenic$ensembl %in% coloc$ensembl))

per_gene_results = coloc %>% group_by(hgnc, ensembl) %>% summarize(best_clpp_mod = max(clpp_mod))

# How many of all tested genes did / didn't colocalize
sum(per_gene_results$best_clpp_mod > 0.4)
sum(per_gene_results$best_clpp_mod <= 0.4)

# How many of monogenic tested genes did / didn't colocalize
monogenic_results = per_gene_results[per_gene_results$ensembl %in% monogenic$ensembl,]
sum(monogenic_results$best_clpp_mod > 0.4)
sum(monogenic_results$best_clpp_mod <= 0.4)

# Print table of results for all monogenic risk genes
print(as.data.frame(monogenic_results %>% arrange(-best_clpp_mod)))
