### All code by Brunilda Balliu 12/21/2019
# Updated by Mike Gloudemans 6/18/2020

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)

#source(file = 'scripts/post_coloc/00_MoMeIR_Functions.R')
save_plots=TRUE

color_scheme = c("white", "#FCD9DD","#F04257", "#9DE1FB", "#26BCF7", "#D690DF", "#BD52CB")
#color_scheme1 = c("white", "#F97B62","#C22506", "#8BE49A", "#196B26", "#CC9A8F", "#4A2922")

# Function to create tileplot of coloc results
plot_coloc_results_function=function(data){
  plot=ggplot(data =  data,
              mapping = aes(x=trait_tissue, y = factor(locus_snp_gene, levels = unique(locus_snp_gene)[length(unique(locus_snp_gene)):1]))) + 
    geom_tile(mapping = aes(fill=coloc_class),color="black") + 
    scale_fill_manual(name="CLPP", values = color_scheme ,drop=FALSE) + 
    geom_text(mapping = aes(label=cross), size=7) +
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
          axis.text.y = element_text(size=12), 
          legend.position = 'top', 
          legend.box = "vertical", 
          legend.text = element_text(size=15), 
          legend.title = element_text(size = 12)) +   
    guides(color=guide_legend(title = "", nrow = 2)) + 
    scale_x_discrete(drop=FALSE)
  return(plot)
}



# Read coloc results
coloc_res=fread(file = 'output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt', sep = '\t', header = T, check.names = F)

coloc_res$hgnc[is.na(coloc_res$hgnc)]=coloc_res$ensembl[is.na(coloc_res$hgnc)]
coloc_res$hgnc[coloc_res$hgnc==""]=coloc_res$ensembl[coloc_res$hgnc==""]

coloc_res$qtl_type = ifelse(grepl("eQTL", coloc_res$eqtl_file), "eqtl", "sqtl")

# Extract the tissue for each entry, regardless of QTL type
coloc_res$tissue=gsub(pattern = "data/eqtls/gtex_v8/",replacement = "", x = coloc_res$eqtl_file)
coloc_res$tissue=gsub(pattern = "data/sqtls/gtex_v8/",replacement = "", x = coloc_res$tissue)
coloc_res$tissue=gsub(pattern = ".allpairs.txt.gz.eQTLs.txt.gz",replacement = "", x = coloc_res$tissue)
coloc_res$tissue=gsub(pattern = ".sQTLs.txt.gz",replacement = "", x = coloc_res$tissue)

# Label each locus with a QTL type, to be used for determining colors later
classes = sapply(1:dim(coloc_res)[1],
      function(x)
      {
	      tmp = coloc_res[x]
	      matches = coloc_res[(coloc_res$locus == tmp$locus) & (coloc_res$tissue == tmp$tissue) & (coloc_res$gwas_trait == tmp$gwas_trait) & (coloc_res$ensembl == tmp$ensembl),]

	      # There are quite a few cases to test for
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.4)) > 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.4)) > 0))
	      {
		      return("strong-both")
	      }
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.4)) > 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.4)) == 0))
	      {
		      return("strong-eqtl")
	      }
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.4)) == 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.4)) > 0))
	      {
		      return("strong-sqtl")
	      }
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.25)) > 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.25)) > 0))
	      {
		      return("weak-both")
	      }
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.25)) == 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.25)) > 0))
	      {
		      return("weak-sqtl")
	      }
	      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= 0.25)) > 0) &&
		      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= 0.25)) == 0))
	      {
		      return("weak-eqtl")
	      }
	      return("none")
      }
)

coloc_res$coloc_class = classes
coloc_res$coloc_class=factor(x = coloc_res$coloc_class, 
levels=c("none", 
         "weak-sqtl",  
         "strong-sqtl",  
         "weak-eqtl",
	 "strong-eqtl",
	 "weak-both",
	 "strong-both"))

coloc_res$base_gwas_file=gsub(pattern = "data/gwas/formatted/sumstats/hg38/",replacement = "", x = coloc_res$base_gwas_file)
coloc_res$base_gwas_file=sapply(X = 1:nrow(coloc_res), FUN = function(X) strsplit(x = coloc_res$base_gwas_file[X], split = '/', fixed = T)[[1]][2])
coloc_res$base_gwas_file=gsub(pattern = ".txt.gz",replacement = "", x = coloc_res$base_gwas_file)
coloc_res$base_gwas_file=gsub(pattern = ".txt",replacement = "", x = coloc_res$base_gwas_file)


# Change eQTL tissue legends
coloc_res$tissue=factor(x = coloc_res$tissue, 
levels=c("Adipose_Subcutaneous", 
         "Adipose_Visceral_Omentum",  
         "Muscle_Skeletal",  
         "Liver",
         "Pancreas"),
labels = c("AdpSQ", "AdpV",  "MuSk", "Liv", "Panc"))

# Annotate resuls by significance of GWAS and eQTL hits
coloc_res$cross=" "
coloc_res$cross_col=NA

# Significant in both
coloc_res$cross_col[coloc_res$min_gwas_pval <= 5e-08 & coloc_res$min_eqtl_pval <= 1e-05] = "GWAS P <= 5e-08; eQTL P <= 1e-05"

# "Borderline" significant GWAS
coloc_res$cross[(coloc_res$min_gwas_pval > 5e-08) & (coloc_res$min_gwas_pval <= 5e-05) & (coloc_res$min_eqtl_pval <= 1e-05)] = "X"
coloc_res$cross_col[(coloc_res$min_gwas_pval > 5e-08) & (coloc_res$min_gwas_pval <= 5e-05) & (coloc_res$min_eqtl_pval <= 1e-05)] = "5e-08 < GWAS P <= 5e-05; eQTL P <= 1e-05"

# Insignificant GWAS
coloc_res$cross[coloc_res$min_gwas_pval > 5e-05] = "Z"
coloc_res$cross_col[coloc_res$min_gwas_pval > 5e-05] = "GWAS P > 5e-05"

coloc_res$cross[coloc_res$min_gwas_pval <= 5e-05 & coloc_res$min_eqtl_pval > 1e-05] = "Z"
coloc_res$cross_col[coloc_res$min_gwas_pval <= 5e-05 & coloc_res$min_eqtl_pval > 1e-05] = "GWAS P <= 5e-05; eQTL P >= 1e-05"

# Less than 20 SNPs around locus
coloc_res$cross[coloc_res$n_snps <= 20] = "S"
coloc_res$cross_col[coloc_res$n_snps <= 20] = "N <=20"

coloc_res$cross_col=factor(coloc_res$cross_col, levels = c("GWAS P <= 5e-08; eQTL P <= 1e-05", "5e-08 < GWAS P <= 5e-05; eQTL P <= 1e-05", "GWAS P > 5e-05", "GWAS P <= 5e-05; eQTL P >= 1e-05", "N <=20"))

coloc_res$step1=factor(coloc_res$step1)
coloc_res$step2=factor(coloc_res$step2)
coloc_res$step3=factor(coloc_res$step3)

# Do one where it's stratified by pancreas presence instead
coloc_res$step1_pancreas = factor(paste(coloc_res$step1, coloc_res$pancreas, sep="-"))
coloc_res$step2_pancreas = factor(paste(coloc_res$step2, coloc_res$pancreas, sep="-"))
coloc_res$step3_pancreas = factor(paste(coloc_res$step3, coloc_res$pancreas, sep="-"))

# Find locus-gene pairs with no coloc at all
# I know this isn't efficient but...it's OK for now
coloc_res$blanks = sapply(1:dim(coloc_res)[1],
      function(x)
      {
	      tmp = coloc_res[x]
	      matches = coloc_res[(coloc_res$locus == tmp$locus) & (coloc_res$ensembl == tmp$ensembl),]
	      
	      if (sum(matches$coloc_class!="none") == 0)
	      {
		return("blank")
	      }
	      else
	      {
		return("okay")
	      }
      }
)

# Change GWAS trait legends
coloc_res$gwas_trait=factor(x = coloc_res$gwas_trait, 
   levels = c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz", "FastInsu_adjBMI_MAGIC_Europeans.txt.gz", "FastGlu_MAGIC_Europeans.txt.gz", "T2D_Mahajan_Europeans.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Spracklen_2020.txt.gz", "WHR-adj-BMI_GIANT_2018.txt.gz", "HDL_GLGC_Expanded.txt.gz", "TG_GLGC_Expanded.txt.gz", "BMI_GIANT_2018.txt.gz"),
   labels = c("ISI_Genesis", "ISI_MAGIC", "FastInsu", "FastGlu", "T2D-Mahajan","T2D-Suzuki", "T2D-Xue", "T2D-Spracklen","WHR","HDL","TG","BMI"))

coloc_res_t2d_collapsed = coloc_res
coloc_res_t2d_collapsed$gwas_trait[as.character(coloc_res$gwas_trait) %in% c("T2D_Mahajan_Europeans.txt.gz", "Type-2-Diabetes_Suzuki_2019.txt.gz", "Type-2-Diabetes_Xue_2018.txt.gz", "Type-2-Diabetes_Spracklen_2020.txt.gz")] = "T2D"
coloc_res_t2d_collapsed$gwas_trait=factor(x = coloc_res_t2d_collapsed$gwas_trait, 
   levels = c("MI_adjBMI_European.txt.gz", "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz", "FastInsu_adjBMI_MAGIC_Europeans.txt.gz", "FastGlu_MAGIC_Europeans.txt.gz", "T2D", "WHR-adj-BMI_GIANT_2018.txt.gz", "HDL_GLGC_Expanded.txt.gz", "TG_GLGC_Expanded.txt.gz", "BMI_GIANT_2018.txt.gz"),
   labels = c("ISI_Genesis", "ISI_MAGIC", "FastInsu", "FastGlu", "T2D","WHR","HDL","TG","BMI"))

#  Sort by rsid, gene, trait, and tissue name
coloc_res=coloc_res[with(coloc_res, order(ensembl, gwas_trait, eqtl_file)),]
coloc_res_t2d_collapsed=coloc_res_t2d_collapsed[with(coloc_res_t2d_collapsed, order(ensembl, gwas_trait, eqtl_file)),]

# Remove instances with multiple SNPs at the same locus
coloc_res=coloc_res[with(coloc_res, order(-clpp_mod)),]
coloc_res=coloc_res[!duplicated(coloc_res[,c("locus", "ensembl", "gwas_trait", "eqtl_file")]),]
coloc_res_t2d_collapsed=coloc_res_t2d_collapsed[with(coloc_res_t2d_collapsed, order(-clpp_mod)),]
coloc_res_t2d_collapsed=coloc_res_t2d_collapsed[!duplicated(coloc_res_t2d_collapsed[,c("locus", "ensembl", "gwas_trait", "eqtl_file")]),]

coloc_res_no_blanks = coloc_res[coloc_res$blanks != "blank",]
coloc_res_t2d_no_blanks = coloc_res_t2d_collapsed[coloc_res_t2d_collapsed$blanks != "blank",]

### Plot coloc results

# Step 1: Plot them, separated by number of candidate genes and number of colocs

plot_out_dir = "output/post_coloc/plots/heatmaps/"

chunk_size = 100
row_height=0.2

plot_heatmap = function(coloc_res, strat_step, out_sub_folder)
{
	print(strat_step)
	dir.create(paste0(plot_out_dir, "/", out_sub_folder), recursive = TRUE, showWarnings=FALSE)
	for (i in 1:length(levels(coloc_res[[strat_step]])))
	{
		tmp_data=coloc_res %>% mutate(trait_tissue=factor(paste(gwas_trait,tissue, sep = '-'), levels = as.vector(t(outer(levels(gwas_trait), levels(tissue), FUN = "paste", sep="-"))))) %>% 
				      mutate(locus_snp_gene=paste(locus, hgnc, sep = '-')) %>% 
				      filter(coloc_res[[strat_step]] == levels(coloc_res[[strat_step]])[i])
		
		row_count = length(unique(tmp_data$locus))

		max_chunk = floor(row_count / chunk_size + 1)

		for (chunk in 1:max_chunk)
		{
			tmp_chunk = tmp_data[tmp_data$locus %in% unique(tmp_data$locus)[(1+(chunk-1)*chunk_size):min(chunk_size*chunk, row_count)],]

			# Genes per locus
			genes_per_locus = tmp_chunk %>% group_by(locus) %>% summarize(genes_at_locus=length(unique(locus_snp_gene))) %>% arrange(locus)

			num_cols = length(levels(tmp_chunk$trait_tissue))
			num_tissues = length(unique(coloc_res$tissue))
			num_vert_bars = num_cols / num_tissues - 1
			num_rows = length(unique(tmp_chunk$locus_snp_gene))

			num_horz_bars = length(unique(tmp_chunk$locus))
			horz_breaks = cumsum(genes_per_locus$genes_at_locus)

			my.vertical.lines<-data.frame(x=seq(5.5, 5.5+num_tissues*(num_vert_bars-1), by = 5), y = rep(0.5, num_vert_bars), 
				xend=seq(5.5, 5.5+num_tissues*(num_vert_bars-1), by = 5), yend = rep(num_rows + 0.5, num_vert_bars))
			
			my.horizontal.lines<-data.frame(x=rep(0.5, num_horz_bars), y=num_rows-horz_breaks+0.5, 
				xend=rep(num_cols+0.5, num_horz_bars), yend=num_rows - horz_breaks+0.5)
			
			plot=plot_coloc_results_function(data = tmp_chunk %>% filter(locus %in% unique(tmp_data$locus)[1:length(unique(tmp_data$locus))]))
			plot = plot + geom_segment(data=my.vertical.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
			plot = plot + geom_segment(data=my.horizontal.lines, aes(x,y,xend=xend, yend=yend), size=0.25, inherit.aes=F)

			if(save_plots) 
			{
				ggsave(filename = paste0(plot_out_dir, '/', out_sub_folder, '/CLPP_group_',levels(coloc_res[[strat_step]])[i],'.part', chunk, '.pdf'), plot = plot, width = 15, height = 5+(row_height*num_rows), limitsize = F)
			}
			
			plot
		}
	}
}

strat_steps = c("step1", "step2", "step3", "step1_pancreas", "step2_pancreas", "step3_pancreas")
out_folders = c("standard", "standard", "standard", "pancreas_stratified", "pancreas_stratified", "pancreas_stratified")
noblanks_out_folders = c("concise", "concise", "concise", "concise_pancreas_stratified", "concise_pancreas_stratified", "concise_pancreas_stratified")

for (i in 1:length(strat_steps))
{
	plot_heatmap(coloc_res, strat_steps[i], out_folders[i])
	plot_heatmap(coloc_res_no_blanks, strat_steps[i], noblanks_out_folders[i])
	plot_heatmap(coloc_res_t2d_collapsed, strat_steps[i], t2d_out_folders[i])
	plot_heatmap(coloc_res_t2d_no_blanks, strat_steps[i], t2d_no_blanks_out_folders[i])
}

# Other things:
# collapse across tissues
# collapse across traits
# collapse across loci
# collapse across tissues and traits
# collapse across tissues and traits and loci

# Since these are all just single steps, I'm sure there's a more generalizable way of doing it...


