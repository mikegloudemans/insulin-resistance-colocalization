import subprocess

trait_map = {"BMI_GIANT_2018.txt.gz": "data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz",
            "FastGlu_MAGIC_Europeans.txt.gz": "data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz",        
            "FastInsu_adjBMI_MAGIC_Europeans.txt.gz": "data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz",        
            "HDL_GLGC_Expanded.txt.gz": "data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz",        
            "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz": "data/gwas/formatted/sumstats/hg38/MAGIC_ISI_Model_2_AgeSexBMI.txt/MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz",        
            "MI_adjBMI_European.txt.gz": "data/gwas/formatted/sumstats/hg38/MI_adjBMI_European/MI_adjBMI_European.txt.gz",        
            "T2D_Mahajan_Europeans.txt.gz": "data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz",        
            "TG_GLGC_Expanded.txt.gz": "data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz",        
            "Type-2-Diabetes_Spracklen_2020.txt.gz": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz",        
            "Type-2-Diabetes_Suzuki_2019.txt.gz": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz",        
            "Type-2-Diabetes_Xue_2018.txt.gz": "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz",        
            "WHR-adj-BMI_GIANT_2018.txt.gz": "data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz"
        }


bests = {}

with open("data/curated_gene_sets/single_coloc_genes_updated.txt") as f:
    for line in f:
        bests[line.strip()] = {}

with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt") as f:
    f.readline()
    for line in f:
        data = line.strip().split()
        hgnc = data[12]
        gwas = data[8]
        clpp_mod = float(data[14])
        if hgnc in bests:
            if gwas not in bests[hgnc]:
                bests[hgnc][gwas] = (data[1], int(data[2]), 0)
            if clpp_mod > bests[hgnc][gwas][2]:
                bests[hgnc][gwas] = (data[1], int(data[2]), clpp_mod)

print bests

with open("effect_sizes.txt", "w") as w:
    w.write("gene\tgwas\tabs_beta\tclpp_mod\n")
    for b in bests:
        for gwas in bests[b]:

            beta_column = subprocess.check_output("zcat {0} | head -n 1".format(trait_map[gwas]), shell=True).strip().split().index("beta")

            out = subprocess.check_output("tabix {0} {1}:{2}-{3}".format(trait_map[gwas], bests[b][gwas][0], bests[b][gwas][1], bests[b][gwas][1]), shell=True).strip().split()
            print b, gwas, out, out[beta_column]

            w.write("{0}\t{1}\t{2}\t{3}\n".format(b, gwas, abs(float(out[beta_column])), bests[b][gwas][2]))
