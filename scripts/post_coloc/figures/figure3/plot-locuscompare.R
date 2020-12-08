require(ggplot2)
require(cowplot)
require(locuscomparer)

# Some added stuff that allows me to modify axes in locuscomparer
source("scripts/post_coloc/figures/figure3/lc_mod.R")

# All of these are T2D GWAS

# These ones are for PLEKHA1
plekha1_sqtl = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Adipose_Subcutaneous_sQTLs_txt_gz/10.122427031.122429624.clu_eqtl_lc_data.txt"
plekha1_sqtl_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Adipose_Subcutaneous_sQTLs_txt_gz/10.122427031.122429624.clu_gwas_lc_data.txt"
plekha1_eqtl = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000107679.14_eqtl_lc_data.txt"
plekha1_eqtl_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000107679.14_gwas_lc_data.txt"
plekha1_eqtl_pancreas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Pancreas_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000107679.14_eqtl_lc_data.txt"
plekha1_eqtl_pancreas_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/T2D_Mahajan_Europeans_txt_gz/T2D_Mahajan_Europeans.txt.gz/10_122433665/Pancreas_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000107679.14_gwas_lc_data.txt"

# These ones are for the BDNF-AS locus
bdnf_as_sqtl = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_sQTLs_txt_gz/11.27640005.27659171.clu_eqtl_lc_data.txt"
bdnf_as_sqtl_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_sQTLs_txt_gz/11.27640005.27659171.clu_gwas_lc_data.txt"
bdnf_as_eqtl = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000245573.7_eqtl_lc_data.txt"
bdnf_as_eqtl_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000245573.7_gwas_lc_data.txt"
lin7c_eqtl = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000148943.11_eqtl_lc_data.txt"
lin7c_eqtl_gwas = "output/post_coloc/plots/locus_compare/finemapping_tmp/2020-12-02_12-22-21.345065/locuscompare/Type-2-Diabetes_Spracklen_2020_txt_gz/Type-2-Diabetes_Spracklen_2020.txt.gz/11_27703198/Adipose_Subcutaneous_allpairs_txt_gz_eQTLs_txt_gz/ENSG00000148943.11_gwas_lc_data.txt"

#######################################################

plot_combo = function(gwas_file, qtl_file, xmax, ymax, title, title2, out_file, use_qtl_snp=TRUE, snp_override="")
{

	gwas_table = read.table(gwas_file)
	snp = gwas_table$V1[which.min(as.numeric(as.character(gwas_table$V2), na.rm=TRUE))]

	if (use_qtl_snp == TRUE)
	{
		qtl_table = read.table(qtl_file)
		snp = qtl_table$V1[which.min(as.numeric(as.character(qtl_table$V2), na.rm=TRUE))]
	}

	if (snp_override != "")
	{
		snp = snp_override
	}

	print(snp)
	
	m = locuscompare_custom(in_fn1 = qtl_file, in_fn2 = gwas_file, xmax=xmax, ymax=ymax, title=title, title2=title2, snp=snp, genome="hg38")
	#m = locuscompare_custom(in_fn1 = edqtl, in_fn2 = edqtl_gwas, xmax=8, ymax=8, title="Spleen edQTL, TNFRSF14", title2="IBD GWAS, Liu et al. 2015", snp=snp, genome="hg38")

	m

	ggsave(out_file, width=12, height=6)
	#ggsave("../output/colocalization/custom_coloc_tests/figure-4b/4b_edQTL.pdf", width=12, height=6)

}

plot_combo(plekha1_sqtl, plekha1_sqtl_gwas, 14, 14, "Spleen eQTL, TYK2", "Lupus GWAS, Bentham et al. 2015", "output/post_coloc/plots/figure3/plekha1_sqtl.pdf", snp_override = "rs2280141")
plot_combo(plekha1_eqtl, plekha1_eqtl_gwas, 14, 14, "Spleen edQTL, TYK2", "Lupus GWAS, Bentham et al. 2015", "output/post_coloc/plots/figure3/plekha1_eqtl.pdf", snp_override = "rs2280141")
plot_combo(plekha1_eqtl_pancreas, plekha1_eqtl_pancreas_gwas, 14, 14, "Spleen sQTL, TYK2", "Lupus GWAS, Bentham et al. 2015", "output/post_coloc/plots/figure3/plekha1_eqtl_pancreas.pdf", snp_override = "rs2280141")

plot_combo(bdnf_as_sqtl, bdnf_as_sqtl_gwas, 17, 17, "Spleen eQTL, TYK2", "PBC GWAS, Cordell et al. 2015", "output/post_coloc/plots/figure3/bdnf_as_sqtl.pdf", snp_override = "rs988748")
plot_combo(bdnf_as_eqtl, bdnf_as_eqtl_gwas, 17, 17, "Spleen edQTL, TYK2", "PBC GWAS, Cordell et al. 2015", "output/post_coloc/plots/figure3/bdnf_as_eqtl.pdf", snp_override = "rs988748")
plot_combo(lin7c_eqtl, lin7c_eqtl_gwas, 17, 17, "Spleen sQTL, TYK2", "PBC GWAS, Cordell et al. 2015", "output/post_coloc/plots/figure3/lin7c_eqtl.pdf", snp_override = "rs988748")
