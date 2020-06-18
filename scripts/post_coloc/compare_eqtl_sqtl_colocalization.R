# We're keeping this simple for now...



# Load in the data
data = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_classification_2020-05-11.txt", sep="\t", header=TRUE, fill=TRUE)

# Of all loci that have both a strong sQTL and a strong eQTL colocalization...
both_coloc = data[(nchar(as.character(data$strong_coloc_eqtl_genes)) > 0) & (nchar(as.character(data$strong_coloc_sqtl_genes)) > 0),]

# How often is there an overlapping gene?
sum(nchar(as.character(both_coloc$strong_coloc_both)) > 0)

# How often is there no overlapping gene?
sum(nchar(as.character(both_coloc$strong_coloc_both)) == 0)

