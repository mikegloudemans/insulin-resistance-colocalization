# Make a config file for running the coloc pipeline (with plot_only mode to get locuscompare plots)

gwas_index = {"Fasting-Glucose": "data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz",
              "Fasting-Insulin-BMI-Adjusted": "data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz",
              "BMI": "data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz",
              "High-Density-Lipoprotein": "data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz",
              "Insulin-Sensitivity-Index-Model-2": "data/gwas/formatted/sumstats/hg38/MAGIC_ISI_Model_2_AgeSexBMI.txt/MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz",
              "M-I": "data/gwas/formatted/sumstats/hg38/MI_adjBMI_European/MI_adjBMI_European.txt.gz",
              "Triglycerides": "data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz",
              "Type-2-Diabetes_Mahajan": "data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz",
              "Type-2-Diabetes_Spracklen": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz",
              "Type-2-Diabetes_Suzuki": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz",
              "Type-2-Diabetes_Xue": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz",
              "Waist-Hip-Ratio-BMI-Adjusted": "data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz"}

loci_to_keep = []
with open("scripts/post_coloc/heatmaps/all-figure-3-loci.txt") as f:
    for line in f:
        loci_to_keep.append(line.strip())

with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt") as f:
    with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/locus_compare_list.txt", "w") as w:
        w.write("chr\tsnp_pos\tsource_file\tlookup_file\tsource_trait\tsource_pvalue\tlookup_pvalue\tlookup_trait\n")
        for line in f:
            data = line.strip().split()

            locus = data[3]
            if locus not in loci_to_keep:
                continue

            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(\
                    data[1],
                    data[2],
                    gwas_index[data[6]],
                    data[9],
                    gwas_index[data[6]],
                    1,
                    1,
                    data[10]))
