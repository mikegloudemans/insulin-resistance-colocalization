require(dplyr)
require(reshape)
require(ggplot2)

single_coloc_genes = read.table("data/curated_gene_sets/single_coloc_genes_updated.txt", header=FALSE)
# Get predicted directionality of effects on each trait
data = read.table("output/post_coloc/de_genes/coloc_directionality.txt", sep="\t", header=TRUE)
data$higher_expression_higher_risk = data$higher_expression_higher_risk == "True"
directions = data %>% filter(clpp_mod > 0.35) %>% filter(hgnc %in% single_coloc_genes[,1]) %>% group_by(hgnc, gwas_file) %>% summarize(hehr = sum(!higher_expression_higher_risk) == 0, helr = sum(higher_expression_higher_risk) == 0)
directions = directions %>% filter(!grepl("ISI", gwas_file)) # Because they're aren't any of these in
sum(!directions$hehr & !directions$helr)

directions$direction = "unclear"
directions$direction[directions$hehr == TRUE] = "higher expression -> higher trait level/risk"
directions$direction[directions$helr == TRUE] = "higher expression -> lower trait level/risk"

#collapse T2D
t2d_studies =c("data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz", 
	       "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz", 
	       "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz",
	       "data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz")

directions$gwas_file = as.character(directions$gwas_file)
directions$gwas_file[directions$gwas_file %in% t2d_studies] = "T2D"
# UBE3C just taking majority vote, all others are unanimous
directions = directions %>% filter(!((hgnc == "UBE3C") & (hehr)))

directions$gwas_file = factor(directions$gwas_file, levels = c(
							       "data/gwas/formatted/sumstats/hg38/MI_adjBMI_European/MI_adjBMI_European.txt.gz", 
							       "data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz", 
							       "data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz",
							       "T2D", 
							       "data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz", 
							       "data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz", 
							       "data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz",
							       "data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz"
							       ),
			      	labels = c("MI", "FastInsu", "FastGlu", "T2D", "WHR", "TG", "HDL", "BMI"))

types = directions %>% group_by(hgnc) %>% summarize(t2d_up = sum(gwas_file %in% c("MI", "ISI", "FastInsu", "FastGlu", "T2D") & hehr) > 0,
						    t2d_down = sum(gwas_file %in% c("MI", "ISI", "FastInsu", "FastGlu", "T2D") & helr) > 0,
					   	    whr_up = sum(gwas_file == "WHR" & hehr) > 0,
					   	    whr_down = sum(gwas_file == "WHR" & helr) > 0,
						    bmi_up = sum(gwas_file == "BMI" & hehr) > 0,
						    bmi_down = sum(gwas_file == "BMI" & helr) > 0,
					   	    tg_hdl_up = sum(xor(hehr & gwas_file == "TG", helr & gwas_file == "HDL")) > 0,
					   	    tg_hdl_down = sum(xor(helr & gwas_file == "TG", hehr & gwas_file == "HDL")) > 0) 

types = types %>% arrange(t2d_up, t2d_down, whr_up, whr_down, tg_hdl_up, tg_hdl_down, bmi_up, bmi_down) %>% # Sort by disease of colocalization
		filter(t2d_up | t2d_down | whr_up | whr_down | tg_hdl_up | tg_hdl_down)	# Don't keep it if only BMI

directions$hgnc = factor(directions$hgnc, levels = types$hgnc)

t2d_only = types %>% filter(t2d_up | t2d_down) %>% select(hgnc)
directions_t2d = directions %>% filter(hgnc %in% t2d_only$hgnc)
whr_only = types %>% filter(!(t2d_up | t2d_down) & (whr_up | whr_down)) %>% select(hgnc)
directions_whr = directions %>% filter(hgnc %in% whr_only$hgnc)
tg_hdl_only = types %>% filter(!(t2d_up | t2d_down | whr_up | whr_down) & (tg_hdl_up | tg_hdl_down))
directions_tg_hdl = directions %>% filter(hgnc %in% tg_hdl_only$hgnc)

plot_figure = function(data) {
	g=ggplot(data = data,
		mapping = aes(y = hgnc, x = gwas_file)) + 
		geom_tile(mapping = aes(fill=direction),color="black") + 
		scale_fill_manual(values = c("#E36414", "#5F0F40", "grey"), drop=FALSE) + 
		theme(axis.title = element_blank(), 
			axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
			axis.text.y = element_text(size=12), 
			legend.position = 'top', 
			legend.box = "vertical", 
			legend.text = element_text(size=5), 
			legend.title = element_text(size = 7),
			) +
		scale_x_discrete(drop=drop)

}
plot_figure(directions)	
ggsave("output/post_coloc/plots/figure4/tiled_coloc_directionality_all.pdf", height=20, width=3.5)
plot_figure(directions_t2d)
ggsave("output/post_coloc/plots/figure4/tiled_coloc_directionality_t2d.pdf", height=6.5, width=2.5)
plot_figure(directions_whr)	
ggsave("output/post_coloc/plots/figure4/tiled_coloc_directionality_whr.pdf", height=9, width=3)
plot_figure(directions_tg_hdl)	
ggsave("output/post_coloc/plots/figure4/tiled_coloc_directionality_tg_hdl.pdf", height=8, width=3)
