require(dplyr)

adipose_whr_data = read.table("output/post_coloc/de_genes/perturbation_by_adipose_whr_coloc.txt", header=TRUE)
liver_lipids_data = read.table("output/post_coloc/de_genes/perturbation_by_liver_lipids_coloc.txt", header=TRUE)
muscle_glucose_data = read.table("output/post_coloc/de_genes/perturbation_by_muscle_glucose_coloc.txt", header=TRUE)


# Should be all fat
unique(adipose_whr_data$gwas_file)
unique(adipose_whr_data$eqtl_file)

adipose_effects = adipose_whr_data %>% filter(pert_tissue == "Fat") %>% group_by(pert_tissue, perturbation, gene) %>% summarize(pert_increases_risk=(sum(perturbation_increases_risk=="True") > 0) && (sum(perturbation_increases_risk=="False")==0), pert_decreases_risk = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") == 0), pert_direction_ambiguous = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") > 0))

adipose_counts = adipose_effects %>% group_by(perturbation) %>% summarize(risk_increasing_gene_counts = sum(pert_increases_risk), risk_decreasing_gene_counts = sum(pert_decreases_risk))
adipose_counts$imbalance_pval_nom = pbinom(pmin(adipose_counts$risk_increasing_gene_counts, adipose_counts$risk_decreasing_gene_counts),adipose_counts$risk_increasing_gene_counts+adipose_counts$risk_decreasing_gene_counts,0.5)*2

# Should be all muscle
unique(muscle_glucose_data$gwas_file)
unique(muscle_glucose_data$eqtl_file)

muscle_effects = muscle_glucose_data %>% filter(pert_tissue == "Muscle") %>% group_by(pert_tissue, perturbation, gene) %>% summarize(pert_increases_risk=(sum(perturbation_increases_risk=="True") > 0) && (sum(perturbation_increases_risk=="False")==0), pert_decreases_risk = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") == 0), pert_direction_ambiguous = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") > 0))

muscle_counts = muscle_effects %>% group_by(perturbation) %>% summarize(risk_increasing_gene_counts = sum(pert_increases_risk), risk_decreasing_gene_counts = sum(pert_decreases_risk))
muscle_counts$imbalance_pval_nom = pbinom(pmin(muscle_counts$risk_increasing_gene_counts, muscle_counts$risk_decreasing_gene_counts),muscle_counts$risk_increasing_gene_counts+muscle_counts$risk_decreasing_gene_counts,0.5)*2

# Should be all liver
unique(liver_lipids_data$gwas_file)
unique(liver_lipids_data$eqtl_file)

liver_effects = liver_lipids_data %>% filter(pert_tissue == "Muscle") %>% group_by(pert_tissue, perturbation, gene) %>% summarize(pert_increases_risk=(sum(perturbation_increases_risk=="True") > 0) && (sum(perturbation_increases_risk=="False")==0), pert_decreases_risk = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") == 0), pert_direction_ambiguous = (sum(perturbation_increases_risk=="False") > 0) && (sum(perturbation_increases_risk=="True") > 0))

liver_counts = liver_effects %>% group_by(perturbation) %>% summarize(risk_increasing_gene_counts = sum(pert_increases_risk), risk_decreasing_gene_counts = sum(pert_decreases_risk))
liver_counts$imbalance_pval_nom = pbinom(pmin(liver_counts$risk_increasing_gene_counts, liver_counts$risk_decreasing_gene_counts),liver_counts$risk_increasing_gene_counts+liver_counts$risk_decreasing_gene_counts,0.5)*2


write.table(liver_counts, "output/post_coloc/de_genes/liver_lipids_pert_direction_counts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(muscle_counts, "output/post_coloc/de_genes/muscle_glucose_pert_direction_counts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(adipose_counts, "output/post_coloc/de_genes/adipose_whr_pert_direction_counts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
