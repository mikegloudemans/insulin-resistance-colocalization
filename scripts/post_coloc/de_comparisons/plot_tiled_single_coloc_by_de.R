require(dplyr)
require(ggplot2)

# TODO: This is screwy, try it again, I'm not using the right column

data = read.table("output/post_coloc/de_genes/perturbation_by_single_coloc.txt", header=TRUE)
data$tissue_by_pert = paste(data$pert_tissue, data$perturbation, sep="_")

key_stats = data[c("tissue_by_pert", "hgnc", "pert_direction")]

key_stats$hgnc = factor(key_stats$hgnc, levels = rev(levels(key_stats$hgnc)))
key_stats$pert_direction = factor(key_stats$pert_direction, levels = c("False", "True"), labels=c("decreases", "increases"))

g=ggplot(data = key_stats,
        mapping = aes(y = hgnc, x = tissue_by_pert)) + 
	geom_tile(mapping = aes(fill=pert_direction),color="black") + 
    	scale_fill_manual(name="perturbation effect on expression", values = c("red", "blue"), drop=FALSE) + 
    	theme(axis.title = element_blank(), 
          	axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
          	axis.text.y = element_text(size=12), 
		legend.position = 'top', 
		legend.box = "vertical", 
		legend.text = element_text(size=15), 
		legend.title = element_text(size = 12)) 
    	#guides(color=guide_legend(title = "", nrow = 2)) + 
    	#scale_x_discrete(drop=FALSE)

ggsave("test.pdf", height=14, width=14)

