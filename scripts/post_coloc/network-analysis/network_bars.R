require(dplyr)
require(reshape2)
require(ggplot2)

secondary = read.table("data/network_analysis/tables/ppi_coloc_diabetes_secondary/ppi_coloc_diabetes_secondary.xlsx_Sheet_1.txt", header=TRUE, sep="\t")
degree_count = secondary %>% group_by(Diabetes_Symbol) %>% summarize(first_degree = length(unique(Coloc_Symbol[Shortest_Path==2])), second_degree = length(unique(Coloc_Symbol[Shortest_Path==3])), third_degree = length(unique(Coloc_Symbol[Shortest_Path==4])), top_three_degrees=first_degree+second_degree+third_degree, top_two_degrees=first_degree+second_degree) %>% arrange(-top_two_degrees)

write.table(degree_count, file="output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/network_analysis/ppi_connectivity_counts.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

degree_count = degree_count[1:10,]
degree_count$Diabetes_Symbol = factor(degree_count$Diabetes_Symbol, levels = rev(degree_count$Diabetes_Symbol))

melted = melt(degree_count)
head(melted)
melted = melt(degree_count[,1:4])
head(melted)
unique(melted$variable
)

g = ggplot(melted, aes(x=Diabetes_Symbol, y=value, fill=variable)) + 
	  geom_col(width=0.75, position = position_dodge(width = 0.75)) +
	  coord_flip()
ggsave("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/network_analysis/diabetes_ppi_counts.pdf", width=6, height=6)

