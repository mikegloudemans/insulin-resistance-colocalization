### All code by Brunilda Balliu 12/21/2019
# Updated by Mike Gloudemans 6/18/2020

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)
require(rjson)

################################
# Key variables
################################

# Files

# Adjustable features
save_plots=TRUE

# Load config info
config_file = commandArgs(trailing=TRUE)[1]
config = fromJSON(file=config_file)
coloc_file = config$coloc_file
plot_out_dir = config$out_folder_prefix

chunk_size = 100
if ("chunk_size" in names(config))
{
	chunk_size = as.numeric(config$chunk_size)
}
row_height=0.2
col_width=0.2

# Load color scheme
source("scripts/post_coloc/figures/color_scheme.R")

################################

main = function()
{
	coloc_res = get_coloc_results(coloc_file)

	coloc_res = coloc_res %>% arrange(-clpp_mod)
	
	# Make an individual split for every stratification wanted.
       	# Specify this in the config file	
	for (strat in config$file_strata)
	{
		coloc_res_tmp = coloc_res
		coloc_res_tmp$split_column = ""
		for (column in strat$split_factors)	
		{
			if (sum(coloc_res_tmp$split_column != "") == 0)
			{
				coloc_res_tmp$split_column = coloc_res_tmp[[column]]
			}
			else
			{
				coloc_res_tmp$split_column = paste(coloc_res_tmp$split_column, coloc_res_tmp[[column]], sep="-")
			}
		}
		coloc_res_tmp$split_column = factor(coloc_res_tmp$split_column)

		# Do collapsing of certain GWAS into a single trait
	        # if desired

		if ("gwas_remapping" %in% names(strat))
		{
			coloc_res_tmp$gwas_label = as.character(coloc_res_tmp$gwas_label)
			for (gwas in names(strat$gwas_remapping$remaps))
			{
				coloc_res_tmp$gwas_label[coloc_res_tmp$gwas_label == gwas] = strat$gwas_remapping$remaps[[gwas]]
			}
			coloc_res_tmp$gwas_label = factor(x=coloc_res_tmp$gwas_label, levels=strat$gwas_remapping$new_order)
		}

		coloc_res_tmp = collapse_axis_factors(coloc_res_tmp)

		# Classify the coloc results to determine what color they'll be in the plot
		coloc_res_tmp$coloc_class = get_heatmap_classes(coloc_res_tmp)

		# Find locus-gene pairs with no coloc at all
		# Label rows as "blank" if they have no matches, so we can leave them out of plots
		coloc_res_tmp$blanks = ""
		for (y_factor in unique(coloc_res_tmp$y_factor))
		{
			locus_matches = coloc_res_tmp[coloc_res_tmp$y_factor == y_factor,]
			if (sum(locus_matches$coloc_class!="none") == 0)
			{
				coloc_res_tmp$blanks[(coloc_res_tmp$y_factor == y_factor)] = "blank"
			}
			else
			{
				coloc_res_tmp$blanks[(coloc_res_tmp$y_factor == y_factor)] = "okay"
			}	
		}

		# Remove rows (gene-locus pairs) with no colocs at all, if desired
		if (strat$concise == "True")
		{
			coloc_res_tmp = coloc_res_tmp %>% filter(blanks != "blank")
		}

		# If there are multiple SNPs in the same plot cell, then
		# just pick the one with the highest CLPP mod
		coloc_res_tmp = coloc_res_tmp %>% arrange(-clpp_mod)
		coloc_res_tmp = coloc_res_tmp[!duplicated(coloc_res_tmp[,c("x_factor", "y_factor")]),]

		coloc_res_tmp = coloc_res_tmp %>% arrange(y_factor)

		if (config$cluster == "True")
		{
			# Binarize cells into colocalize or non-colocalized
			coloc_res_tmp$clust_stat = 0
			coloc_res_tmp$clust_stat[coloc_res_tmp$coloc_class == "none"] = 1

			dc = dcast(coloc_res_tmp, y_factor ~ x_factor, value.var = "clust_stat")
			grid = as.matrix(dc[,-1])
			rownames(grid) = dc[,1]
			grid[is.na(grid)] = 0
			
			dst = dist(dc, method = "binary")
			h = hclust(dst)
			
			coloc_res_tmp$y_factor = factor(coloc_res_tmp$y_factor, levels = rownames(grid)[h$order])
		}

		### Plot coloc results
		plot_heatmap(coloc_res_tmp, strat$out_dir)
	}

}

######################################################

# Function to create tileplot of coloc results
plot_coloc_results_function=function(data){
  plot=ggplot(data =  data,
              mapping = aes(x=x_factor, y = y_factor)) + 
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

get_coloc_results = function(coloc_file)
{
	# Read coloc results
	coloc_res=as.data.frame(fread(file = coloc_file, sep = '\t', header = T, check.names = F))

	# NOTE: these lines are no longer necessary if the latest version of aggregate_coloc_results.R
	# was used to generate the input file
	#coloc_res$hgnc[is.na(coloc_res$hgnc)]=coloc_res$ensembl[is.na(coloc_res$hgnc)]
	#coloc_res$hgnc[coloc_res$hgnc==""]=coloc_res$ensembl[coloc_res$hgnc==""]

	# Identify the QTL type and tissue for each coloc test
	coloc_res$qtl_type = ""
	coloc_res$tissue = ""
	for (qtl_type in names(config$qtl_types))
	{
		for (file in names(config$qtl_types[[qtl_type]]))
		{
			coloc_res$qtl_type[coloc_res$eqtl_file == file] = qtl_type
			coloc_res$tissue[coloc_res$eqtl_file == file] = config$qtl_types[[qtl_type]][[file]]
		}
	}
	coloc_res$tissue = factor(x=coloc_res$tissue, levels=config$tissue_order)

	# Label the GWAS'es with short and readable tags
	coloc_res$gwas_label = ""
	for (file in names(config$gwas_labels))
	{
		coloc_res$gwas_label[coloc_res$base_gwas_file == file] = config$gwas_labels[[file]]
	}
	coloc_res$gwas_label = factor(x=coloc_res$gwas_label, levels=config$gwas_order)

	if (("mark_dubious_results" %in% names(config)) && (config$mark_dubious_results == "True"))
	{
		# IF annotating with crosses and other marks to denote which ones 
		# aren't trustworthy, do it here.
		coloc_res = mark_dubious_results(coloc_res)
	}
	else
	{
		coloc_res$cross = ""
		coloc_res$cross_col = ""
		coloc_res$cross = factor(coloc_res$cross)
		coloc_res$cross_col = factor(coloc_res$cross_col)
	}

	return(coloc_res)
}

mark_dubious_results = function(coloc_res)
{
	# TODO: Make this more general

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
} 

plot_heatmap = function(coloc_res, out_sub_folder)
{
	dir.create(paste0(plot_out_dir, "/", out_sub_folder), recursive = TRUE, showWarnings=FALSE)
	splits = unique(levels(coloc_res[["split_column"]]))
	for (i in 1:length(splits))
	{
		split_col = splits[i]
		print(paste0("Plotting by ", split_col))

		tmp_data=coloc_res %>% 
			filter(coloc_res[["split_column"]] == splits[i])
		
		row_count = length(unique(tmp_data$y_factor))

		max_chunk = floor(row_count / chunk_size + 1)

		for (chunk in 1:max_chunk)
		{
			tmp_chunk = tmp_data[tmp_data$y_factor %in% unique(tmp_data$y_factor)[(1+(chunk-1)*chunk_size):min(chunk_size*chunk, row_count)],]

			# Genes per locus
			genes_per_locus = tmp_chunk %>% group_by(locus) %>% summarize(genes_at_locus=length(unique(y_factor))) %>% arrange(locus)

			num_cols = length(levels(tmp_chunk$x_factor))
			if (("x_axis_collapse" %in% names(config)) && ((config$x_axis_collapse == "tissues") || (config$x_axis_collapse == "tissues-gwas")))
			{
				num_tissues = 1
			} else
			{
				num_tissues = length(unique(coloc_res$tissue))
			}
			num_vert_bars = num_cols / num_tissues - 1
			num_rows = length(unique(tmp_chunk$y_factor))

			num_horz_bars = length(unique(tmp_chunk$locus))
			horz_breaks = cumsum(genes_per_locus$genes_at_locus)

			margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$y_factor))), 0))
			
			plot=plot_coloc_results_function(data = tmp_chunk)
			
			if (num_vert_bars != 0)
			{
				my.vertical.lines<-data.frame(x=seq(0, num_vert_bars-1, by = num_tissues) + num_tissues + 0.5, y = rep(0.5, num_vert_bars), 
					xend=seq(0, num_vert_bars-1, by = num_tissues) + num_tissues + 0.5, yend = rep(num_rows + 0.5, num_vert_bars))
				plot = plot + geom_segment(data=my.vertical.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
			}
			if (num_horz_bars != 0)
			{
				my.horizontal.lines<-data.frame(x=rep(0.5, num_horz_bars), y=num_rows-horz_breaks+0.5, 
					xend=rep(num_cols+0.5, num_horz_bars), yend=num_rows - horz_breaks+0.5)
				plot = plot + geom_segment(data=my.horizontal.lines, aes(x,y,xend=xend, yend=yend), size=0.25, inherit.aes=F)
			}
			plot

			if(save_plots) 
			{
				ggsave(filename = paste0(plot_out_dir, '/', out_sub_folder, '/CLPP_group_',levels(coloc_res[["split_column"]])[i],'.part', chunk, '.pdf'), plot = plot, width = margin_approx_size+(col_width*num_cols), height = margin_approx_size+(row_height*num_rows), limitsize = F)
			}
			
		}
	}
}

get_heatmap_classes = function(coloc_res)
{
	# TODO: I need to generalize this to other things beyond eqtl and sqtl
	# Label each locus with a QTL type, to be used for determining colors later
	classes = sapply(1:dim(coloc_res)[1],
	      function(x)
	      {
		      tmp = coloc_res[x,]
		      matches = coloc_res[(coloc_res$x_factor == tmp$x_factor) & (coloc_res$y_factor == tmp$y_factor),]

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
}

collapse_axis_factors = function(coloc_file)
{
	coloc_res = coloc_file

	# Collapse genes by locus
	# Remove all but the best coloc at each locus, regardless of gene
	if (("y_axis_collapse" %in% names(config)) && (config$y_axis_collapse == "genes"))
	{
		coloc_res$y_factor = coloc_res$locus
	} else
	{
		coloc_res$y_factor = paste(coloc_res$locus, coloc_res$hgnc, sep="-")
	}

	coloc_res = coloc_res %>% arrange(as.numeric(coloc_res$locus))
	coloc_res$y_factor = factor(coloc_res$y_factor, levels = rev(unique(coloc_res$y_factor)))
	
	# Collapse across tissues
	# Remove all but the best tissue for a given
	if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "tissues"))
	{
		coloc_res$x_factor = coloc_res$gwas_label
	} else if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "gwas"))
	{
		coloc_res$x_factor = coloc_res$tissue
	} else if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "tissues-gwas"))
	{
		coloc_res$x_factor = "any"
		coloc_res$x_factor = factor(coloc_res$x_factor)
	} else
	{
		coloc_res$x_factor = paste(coloc_res$gwas_label, coloc_res$tissue, sep="-")
		coloc_res$x_factor = factor(coloc_res$x_factor, levels = as.vector(t(outer(levels(coloc_res$gwas_label), levels(coloc_res$tissue), FUN = "paste", sep="-"))))
	}
	# No need to sort on the x-axis, because the desired orders have already been pre-specified
	
	return(coloc_res)
}

main()
