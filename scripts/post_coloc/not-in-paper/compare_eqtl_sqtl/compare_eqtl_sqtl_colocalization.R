# We're keeping this simple for now...

# Load in the data
data = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_classification_2020-05-11.txt", sep="\t", header=TRUE, fill=TRUE)

no_coloc = data[(nchar(as.character(data$strong_coloc_eqtl_genes)) == 0) & (nchar(as.character(data$strong_coloc_sqtl_genes)) == 0),]
eqtl_coloc = data[(nchar(as.character(data$strong_coloc_eqtl_genes)) > 0) & (nchar(as.character(data$strong_coloc_sqtl_genes)) == 0),]
sqtl_coloc = data[(nchar(as.character(data$strong_coloc_eqtl_genes)) == 0) & (nchar(as.character(data$strong_coloc_sqtl_genes)) > 0),]
print("Tested loci with no coloc at all:")
print(dim(no_coloc)[1])
print("Loci with eQTL coloc only:")
print(dim(eqtl_coloc)[1])
print("Loci with sQTL coloc only:")
print(dim(sqtl_coloc)[1])

# Of all loci that have both a strong sQTL and a strong eQTL colocalization...
print("Loci with coloc in both:")
both_coloc = data[(nchar(as.character(data$strong_coloc_eqtl_genes)) > 0) & (nchar(as.character(data$strong_coloc_sqtl_genes)) > 0),]
print(dim(both_coloc)[1])

# How often is there an overlapping gene?
print("Overlapping gene:")
sum(nchar(as.character(both_coloc$strong_coloc_both)) > 0)

# How often is there no overlapping gene?
print("No overlapping gene:")
sum(nchar(as.character(both_coloc$strong_coloc_both)) == 0)

