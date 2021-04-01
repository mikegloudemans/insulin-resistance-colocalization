
# Do each of these steps on a gene-by-gene level, but also in a more aggregated way -- e.g. overall enrichments for all tissue-specific colocs of a certain group or something

gene_set = []
with open("data/curated_gene_sets/single_coloc_genes.txt") as f:
    for line in f:
        gene_set.append(line.strip())


# Plot KO enrichments

# Plot PheWAS enrichments

# Plot chem enrichments


