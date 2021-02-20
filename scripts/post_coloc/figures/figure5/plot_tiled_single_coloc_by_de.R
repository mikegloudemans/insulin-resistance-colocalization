require(dplyr)
require(ggplot2)

# TODO: This is screwy, try it again, I'm not using the right column

data = read.table("output/post_coloc/de_genes/perturbation_by_single_coloc.txt", header=TRUE)
data$tissue_by_pert = paste(data$pert_tissue, data$perturbation, sep="_")

key_stats = data[c("tissue_by_pert", "hgnc", "pert_direction", "pert_tissue", "perturbation")]

key_stats$hgnc = factor(key_stats$hgnc, levels = rev(levels(key_stats$hgnc)))
key_stats$pert_direction = factor(key_stats$pert_direction, levels = c("False", "True"), labels=c("decreases", "increases"))

table(key_stats$pert_tissue, key_stats$perturbation)
colSums(table(key_stats$pert_tissue, key_stats$perturbation))
sort(colSums(table(key_stats$pert_tissue, key_stats$perturbation)))

accepted_perturbations = c("TGFB1",  "RETA", "IGF1",  "IL-6", "SP60",  "GLUC", "TNFa", "IBMX", "INSU", "DEXA")


### Part 1: Tissues ###

tissue_perturbations = c("IL-6", "ISOP", "LAUR", "ROSI")

tissue_stats = key_stats %>% filter(perturbation %in% tissue_perturbations)

gene_order = sapply(unique(tissue_stats$hgnc), function(x)
        {
		sig_tissues = c()
		for (tissue in c("Fat", "Liver", "Muscle"))
		{
			if (sum((tissue_stats$hgnc == x) & (tissue_stats$pert_tissue == tissue)) > sum(tissue_stats$hgnc == x) * 0)
			{
				sig_tissues = c(sig_tissues, tissue)
			}
		}       
		if (length(sig_tissues) == 0)
		{
			return("none")
		} else
		{
			return(paste(sig_tissues, collapse="-"))
		}
	})
names(gene_order) = unique(tissue_stats$hgnc)

print(table(gene_order))
tissue_stats$hgnc = factor(tissue_stats$hgnc, levels = rev(names(sort(gene_order))))

g=ggplot(data = tissue_stats,
        mapping = aes(y = hgnc, x = tissue_by_pert)) + 
	geom_tile(mapping = aes(fill=pert_direction),color="black") + 
    	scale_fill_manual(name="perturbation effect on expression", values = c("#008080", "#FFAC00"), drop=FALSE) + 
    	theme(axis.title = element_blank(), 
          	axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
          	axis.text.y = element_text(size=12), 
		legend.position = 'top', 
		legend.box = "vertical", 
		legend.text = element_text(size=15), 
		legend.title = element_text(size = 12)) 

ggsave("output/post_coloc/plots/figure5/tiled_de_by_coloc_tissues.pdf", height=18, width=6)


### Part 2: Tissues ###

top_perturbations = c("DEXA", "GLUC", "IGF1", "INSU", "TNFa")

top_stats = key_stats %>% filter(perturbation %in% top_perturbations)

gene_order = unique(top_stats$hgnc)
pert_order = sapply(unique(top_stats$hgnc), function(x)
        {
		sig_perts = c()
		for (pert in top_perturbations)
		{
			if (sum((top_stats$hgnc == x) & (top_stats$perturbation == pert)) > sum(top_stats$hgnc == x) * 0)
			{
				sig_perts = c(sig_perts, pert)
			}
		}       
		if (length(sig_perts) == 0)
		{
			return("none")
		} else
		{
			return(paste(sig_perts, collapse="-"))
		}
	})

pert_order = factor(pert_order, levels=rev(c("GLUC",
		       "GLUC-IGF1",
		       "DEXA-GLUC-TNFa",
		       "GLUC-INSU",
		       "GLUC-IGF1-INSU",
		       "GLUC-IGF1-INSU-TNFa",
		       "DEXA-GLUC-INSU-TNFa",
		       "DEXA-GLUC-IGF1-INSU-TNFa",
		       "DEXA-GLUC-IGF1-INSU",
		       "INSU",
		       "DEXA-INSU",
		       "DEXA-INSU-TNFa",
		       "DEXA-IGF1-INSU-TNFa",
		       "DEXA-IGF1-INSU",
		       "IGF1-INSU",
		       "IGF1-INSU-TNFa",
		       "IGF1",
		       "DEXA-IGF1",
		       "TNFa",
		       "DEXA-TNFa",
		       "DEXA"
		       )))

gene_order = gene_order[order(pert_order)]

top_stats$hgnc = factor(top_stats$hgnc, levels = gene_order)

top_stats$tissue_by_pert = factor(top_stats$tissue_by_pert, levels = c("Fat_GLUC", "Liver_GLUC", "Muscle_GLUC", 
					     "Fat_INSU", "Liver_INSU", "Muscle_INSU",
					     "Fat_IGF1", "Liver_IGF1", "Muscle_IGF1",
					     "Fat_TNFa", "Liver_TNFa", "Muscle_TNFa",
					     "Fat_DEXA", "Liver_DEXA","Muscle_DEXA"))

g=ggplot(data = top_stats,
        mapping = aes(y = hgnc, x = tissue_by_pert)) + 
	geom_tile(mapping = aes(fill=pert_direction),color="black") + 
    	scale_fill_manual(name="perturbation effect on expression", values = c("#008080", "#FFAC00"), drop=FALSE) + 
    	theme(axis.title = element_blank(), 
          	axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
          	axis.text.y = element_text(size=12), 
		legend.position = 'top', 
		legend.box = "vertical", 
		legend.text = element_text(size=15), 
		legend.title = element_text(size = 12)) 
    	#guides(color=guide_legend(title = "", nrow = 2)) + 
    	#scale_x_discrete(drop=FALSE)

ggsave("output/post_coloc/plots/figure5/tiled_de_by_coloc_perturbations.pdf", height=18, width=6)

