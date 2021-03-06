### All code by Brunilda Balliu 12/21/2019
# Updated by Mike Gloudemans 6/18/2020

suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(require(rjson)))

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
if ("chunk_size" %in% names(config))
{
	chunk_size = as.numeric(config$chunk_size)
}
row_height=0.21
col_width=0.2

# Load color scheme
source("scripts/post_coloc/figures/color_scheme.R")

################################

main = function()
{	
	# Catch config errors at the beginning rather than after wasting a lot of time
	validate_config(config)

	coloc_res = get_coloc_results(coloc_file)

	coloc_res = coloc_res %>% arrange(-clpp_mod)
	
	# Make an individual split for every stratification wanted.
       	# Specify this in the config file	
	for (strat in config$file_strata)
	{
		coloc_res_tmp = coloc_res
		coloc_res_tmp$split_column = ""
		if(!("split_factors" %in% names(strat)))
		{
			coloc_res_tmp$split_column = "no-split"
		}
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

		if ("gwas_blacklist" %in% names(strat))
		{
			for (bl in strat$gwas_blacklist)
			{
				coloc_res_tmp = coloc_res_tmp[coloc_res_tmp$gwas_label != bl,]
				strat$gwas_remapping$new_order = strat$gwas_remapping$new_order[strat$gwas_remapping$new_order != bl]
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
			
			dst = suppressWarnings(dist(dc, method = "binary"))
			dst[is.na(dst)] = 1
			h = hclust(dst)

			#ixs = sort(h$order, index.return=TRUE)$ix

			coloc_res_tmp$y_factor = factor(coloc_res_tmp$y_factor, levels = rownames(grid)[h$order])
		}

		### Plot coloc results
		plot_heatmap(coloc_res_tmp, strat)
	}

}

######################################################

# Function to create tileplot of coloc results
plot_coloc_results_function=function(data){
  plot=ggplot(data =  data,
              mapping = aes(x=x_factor, y = y_factor)) + 
    geom_tile(mapping = aes(fill=coloc_class),color="black") + 
    scale_fill_manual(name="CLPP", values = color_scheme ,drop=FALSE) + 
    geom_text(mapping = aes(label=cross, color=cross_col), size=4) +
    scale_colour_manual(values=c("black", "white")) +
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

	if (("label_individual_cells" %in% names(config)) && (config$label_individual_cells == "True"))
	{
		# IF annotating with crosses and other marks to denote which ones 
		# aren't trustworthy, do it here.
		coloc_res = label_individual_cells(coloc_res)
	} else if (("put_scores_in_cells" %in% names(config)) && (config$put_scores_in_cells) == "True")
	{
		coloc_res = put_scores_in_cells(coloc_res)
	} else
	{
		coloc_res$cross = ""
		coloc_res$cross_col = ""
	}

	coloc_res$cross = factor(coloc_res$cross)
	coloc_res$cross_col = factor(coloc_res$cross_col)

	return(coloc_res)
}

put_scores_in_cells = function(coloc_res)
{
	# Annotate results by significance of GWAS and eQTL hits

	coloc_res_tmp = coloc_res
	coloc_res_tmp$cross = as.character(round(coloc_res_tmp$clpp_mod, 2) * 100)
	coloc_res_tmp$cross_col = "black"
	# NOTE: in the long term we don't want this to be hard-coded
	coloc_res_tmp$cross_col[coloc_res_tmp$clpp_mod > 0.35] = "white"

	return(coloc_res_tmp)
} 

label_individual_cells = function(coloc_res)
{
	# Annotate results by significance of GWAS and eQTL hits

	coloc_res_tmp = coloc_res

	coloc_res_tmp$cross=" "
	coloc_res_tmp$cross_col=NA

	for (i in 1:length(config$label_cells_key))
	{
		class = config$label_cells_key[[i]]

		subset = rep(TRUE, dim(coloc_res_tmp)[1])
		
		if ("gwas_max" %in% names(class))
		{
			subset[coloc_res_tmp$min_gwas_pval < as.numeric(class$gwas_max)] = FALSE
		}	
		if ("gwas_min" %in% names(class))
		{
			subset[coloc_res_tmp$min_gwas_pval >= as.numeric(class$gwas_min)] = FALSE
		}	
		if ("eqtl_max" %in% names(class))
		{
			subset[coloc_res_tmp$min_eqtl_pval < as.numeric(class$eqtl_max)] = FALSE
		}	
		if ("eqtl_min" %in% names(class))
		{
			subset[coloc_res_tmp$min_eqtl_pval >= as.numeric(class$eqtl_min)] = FALSE
		}	
		if ("snp_max" %in% names(class))
		{
			subset[coloc_res_tmp$n_snps < as.numeric(class$snp_max)] = FALSE
		}	
		if ("snp_min" %in% names(class))
		{
			subset[coloc_res_tmp$n_snps >= as.numeric(class$snp_min)] = FALSE
		}
		coloc_res_tmp$cross[subset] = class$mark_char
		coloc_res_tmp$cross_col[subset] = class$caption
	}

	return(coloc_res_tmp)
} 

plot_heatmap = function(coloc_res, strat)
{
	out_sub_folder = strat$out_dir
	dir.create(paste0(plot_out_dir, "/", out_sub_folder), recursive = TRUE, showWarnings=FALSE)
	if ("constrain_split" %in% names(strat))
	{
		splits = strat$constrain_split
	} else
	{
		splits = unique(levels(coloc_res[["split_column"]]))
	}
	for (i in 1:length(splits))
	{
		split_col = splits[i]
		
		print(paste0("Plotting by ", split_col))

		tmp_data=coloc_res %>% 
			filter(coloc_res[["split_column"]] == split_col)
		print(dim(tmp_data))

		row_count = length(unique(tmp_data$y_factor))

		max_chunk = floor(row_count / chunk_size + 1)

		for (chunk in 1:max_chunk)
		{
			tmp_chunk = tmp_data[tmp_data$y_factor %in% unique(tmp_data$y_factor)[(1+(chunk-1)*chunk_size):min(chunk_size*chunk, row_count)],]
			tmp_chunk$y_factor = factor(tmp_chunk$y_factor, levels = rev(unique(tmp_chunk$y_factor)))

			# Genes per locus
			genes_per_locus = tmp_chunk %>% group_by(locus) %>% summarize(genes_at_locus=length(unique(y_factor))) %>% arrange(locus)

			num_cols = length(levels(tmp_chunk$x_factor))
			if (("x_axis_collapse" %in% names(config)) && ((config$x_axis_collapse == "tissues") || (config$x_axis_collapse == "tissues-gwas")))
			{
				num_tissues = 1
			} else
			{
				num_tissues = length(levels(coloc_res$tissue))
			}
			num_vert_bars = num_cols / num_tissues - 1
			num_rows = length(unique(tmp_chunk$y_factor))

			if ((("cluster" %in% names(config)) && (config$cluster == "True")) ||
			    (("y_axis_collapse" %in% names(config)) && (config$y_axis_collapse == "genes")))
			{
				# It doesn't make sense to separate loci if they're clustered
				# It also doesn't make sense to separate loci if each row is a distinct locus
				num_horz_bars = 0
			} else 
			{
				num_horz_bars = length(unique(tmp_chunk$locus))
			}
			horz_breaks = cumsum(genes_per_locus$genes_at_locus)

			y_margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$y_factor))), 0))
			x_margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$x_factor))), 0))+1.3

			plot=plot_coloc_results_function(data = tmp_chunk)
			
			if (num_vert_bars != 0)
			{
				my.vertical.lines<-data.frame(x=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, y = rep(0.5, num_vert_bars), 
					xend=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, yend = rep(num_rows + 0.5, num_vert_bars))
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
				ggsave(filename = paste0(plot_out_dir, '/', out_sub_folder, '/CLPP_group_',split_col,'.part', chunk, '.pdf'), plot = plot, width = y_margin_approx_size+(col_width*num_cols), height = x_margin_approx_size+(row_height*num_rows), limitsize = F)
			}
			
		}
	}
}

get_heatmap_classes = function(coloc_res)
{
	strong_threshold = 0.35
	weak_threshold = 0.35

	# TODO: I need to generalize this to other things beyond eqtl and sqtl
	# Label each locus with a QTL type, to be used for determining colors later
	classes = sapply(1:dim(coloc_res)[1],
	      function(x)
	      {
		      tmp = coloc_res[x,]
		      matches = coloc_res[(coloc_res$x_factor == tmp$x_factor) & (coloc_res$y_factor == tmp$y_factor),]

		      # There are quite a few cases to test for
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= strong_threshold)) > 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= strong_threshold)) > 0))
		      {
			      return("both")
		      }
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= strong_threshold)) > 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= strong_threshold)) == 0))
		      {
			      return("eqtl")
		      }
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= strong_threshold)) == 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= strong_threshold)) > 0))
		      {
			      return("sqtl")
		      }
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= weak_threshold)) > 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= weak_threshold)) > 0))
		      {
			      return("both")
		      }
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= weak_threshold)) == 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= weak_threshold)) > 0))
		      {
			      return("sqtl")
		      }
		      if ((sum((matches$qtl_type == "eqtl") & (matches$clpp_mod >= weak_threshold)) > 0) &&
			      (sum((matches$qtl_type == "sqtl") & (matches$clpp_mod >= weak_threshold)) == 0))
		      {
			      return("eqtl")
		      }
		      return("none")
	      }
	)

	coloc_res$coloc_class = classes
	
	coloc_res$coloc_class=factor(x = coloc_res$coloc_class, 
	levels=c("none", 
		 "sqtl",  
		 "eqtl",
		 "both"))
}

collapse_axis_factors = function(coloc_file)
{
	coloc_res = coloc_file

	if ("locus_selection_list" %in% names(config))
	{
		coloc_res = coloc_res %>% filter(locus %in% config$locus_selection_list)
	}
	if ("gene_selection_list" %in% names(config))
	{
		coloc_res = coloc_res %>% filter(ensembl %in% config$gene_selection_list)
	}

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
	coloc_res$y_factor = factor(coloc_res$y_factor, levels = unique(coloc_res$y_factor))
	
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

validate_config = function(config)
{

	if ( ("label_individual_cells" %in% names(config)) && 
	    (config$label_individual_cells == "True") && 
	    (!(("label_cells_key") %in% names(config))))
	{
		print("'label_individual_cells' requested in config, but no 'label_cells_key' specified.")
		print("Either specify 'label_cells_key' or set 'label_individual_cells' to 'False'.")
		stopifnot(FALSE)
	}
}

main()
