require(dplyr)

gene_list = read.table("data/curated_gene_sets/single_coloc_genes_updated.txt", header=FALSE)
colnames(gene_list) = "Gene"

coloc_data = read.table("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt", header=TRUE, sep="\t")
de_data = read.table("output/post_coloc/de_genes/perturbation_by_single_coloc.txt", header=TRUE, sep="\t")
network_data = read.table("data/network_analysis/tables/ppi_coloc_diabetes_secondary/ppi_coloc_diabetes_secondary.xlsx_Sheet_1.txt", header=TRUE, sep="\t")

add_priority_scores = function(gene_list, coloc_data, de_data, network_data, suffix)
{
	gene_list_copy = gene_list

	for (i in 1:dim(gene_list)[1])
	{
		this_gene = gene_list[i,1]

		# Get colocalization info

		all_coloc = coloc_data %>% filter(as.character(hgnc) == as.character(this_gene))
		t2d_coloc = all_coloc %>% filter(gwas_short %in% c("Fasting-Glucose", "Fasting-Insulin-BMI-Adjusted", "Insulin-Sensitivity-Index-Model-2", "M-I", "Type-2-Diabetes_Mahajan", "Type-2-Diabetes_Spracklen", "Type-2-Diabetes_Suzuki", "Type-2-Diabetes_Xue"))
		whr_coloc = all_coloc %>% filter(gwas_short %in% c("Waist-Hip-Ratio-BMI-Adjusted"))
		tg_hdl_coloc = all_coloc %>% filter(gwas_short %in% c("Triglycerides", "High-Density-Lipoprotein"))
		bmi_coloc = all_coloc %>% filter(gwas_short %in% c("BMI"))
		pancreas_coloc = all_coloc %>% filter(grepl("Pancreas", eqtl_short))

		gene_list_copy[[paste0("has_t2d_coloc", suffix)]][i] = dim(t2d_coloc)[1] > 0
		gene_list_copy[[paste0("has_whr_coloc", suffix)]][i] = dim(whr_coloc)[1] > 0
		gene_list_copy[[paste0("has_tg_hdl_coloc", suffix)]][i] = dim(tg_hdl_coloc)[1] > 0
		gene_list_copy[[paste0("has_bmi_coloc", suffix)]][i] = dim(bmi_coloc)[1] > 0
		gene_list_copy[[paste0("has_pancreas_coloc", suffix)]][i] = dim(pancreas_coloc)[1] > 0
		gene_list_copy[[paste0("all_coloc_score", suffix)]][i] = ((dim(t2d_coloc)[1] > 0) + (dim(whr_coloc)[1] > 0) + (dim(tg_hdl_coloc)[1] > 0) + (dim(bmi_coloc)[1] > 0)) / ifelse(gene_list_copy[["has_pancreas_coloc_any"]][i],2,1) / 4

		# Get perturbation DE info
		all_de = de_data %>% filter(as.character(hgnc) == as.character(this_gene))
		gene_list_copy[[paste0("de_conditions", suffix)]][i] = dim(all_de)[1]
		gene_list_copy[[paste0("has_de", suffix)]][i] =  gene_list_copy[[paste0("de_conditions", suffix)]][i] > 0
		if (suffix == "_any")
		{
			gene_list_copy[[paste0("has_many_de", suffix)]][i] =  gene_list_copy[[paste0("de_conditions", suffix)]][i] > 9
		} else
		{
			gene_list_copy[[paste0("has_many_de", suffix)]][i] =  gene_list_copy[[paste0("de_conditions", suffix)]][i] > 3
		}
			
		gene_list_copy[[paste0("all_de_score", suffix)]][i] = 0
		if (gene_list_copy[[paste0("has_de", suffix)]][i])
		{
			gene_list_copy[[paste0("all_de_score", suffix)]][i] = 0.5
		} 
		if (gene_list_copy[[paste0("has_many_de", suffix)]][i])
		{
			gene_list_copy[[paste0("all_de_score", suffix)]][i] = 1
		}

		# Get network info
		gene_list_copy[[paste0("network_score", suffix)]][i] = 0
		all_network = network_data %>% filter(as.character(Coloc_Symbol) == as.character(this_gene) & Shortest_Path <= 3)	
		if (dim(all_network)[1] > 0)
		{
			gene_list_copy[[paste0("network_score", suffix)]][i] = 1/3
		}
		if (dim(all_network)[1] > 1)
		{
			gene_list_copy[[paste0("network_score", suffix)]][i] = 2/3
		}
		first_order_all_network = all_network %>% filter(Shortest_Path == 2)
		if (dim(first_order_all_network)[1] > 0)
		{
			gene_list_copy[[paste0("network_score", suffix)]][i] = 1
		}

		gene_list_copy[[paste0("score", suffix)]][i] = gene_list_copy[[paste0("all_coloc_score", suffix)]][i] + gene_list_copy[[paste0("all_de_score", suffix)]][i] + gene_list_copy[[paste0("network_score", suffix)]][i] 
	}

	return(gene_list_copy)
}


adipose_coloc_data = coloc_data %>% filter(grepl("Adipose", eqtl_short))
liver_coloc_data = coloc_data %>% filter(grepl("Liver", eqtl_short))
muscle_coloc_data = coloc_data %>% filter(grepl("Muscle", eqtl_short))

adipose_de_data = de_data %>% filter(pert_tissue == "Fat")
liver_de_data = de_data %>% filter(pert_tissue == "Liver")
muscle_de_data = de_data %>% filter(pert_tissue == "Muscle")

adipose_network_data = network_data %>% filter(grepl("Fat", Perturbation))
liver_network_data = network_data %>% filter(grepl("Liver", Perturbation))
muscle_network_data = network_data %>% filter(grepl("Muscle", Perturbation))

gene_list = add_priority_scores(gene_list, coloc_data, de_data, network_data, "_any")
gene_list = add_priority_scores(gene_list, adipose_coloc_data, adipose_de_data, adipose_network_data, "_adipose")
gene_list = add_priority_scores(gene_list, liver_coloc_data, liver_de_data, liver_network_data, "_liver")
gene_list = add_priority_scores(gene_list, muscle_coloc_data, muscle_de_data, muscle_network_data, "_muscle")

gene_list$tissue = "Shared"
gene_list$tissue[(gene_list$score_adipose - gene_list$score_liver >= 1) & (gene_list$score_adipose - gene_list$score_muscle >= 1)] = "Adipose"
gene_list$tissue[(gene_list$score_liver - gene_list$score_adipose >= 1) & (gene_list$score_liver - gene_list$score_muscle >= 1)] = "Liver"
gene_list$tissue[(gene_list$score_muscle - gene_list$score_liver >= 1) & (gene_list$score_muscle - gene_list$score_adipose >= 1)] = "Muscle"

write.table(gene_list, file = "output/post_coloc/plots/figure6/gene_priority_list.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
