import glob
import subprocess
import gzip

gwas_aliases = \
{
    "WHR-adj-BMI_GIANT_2018_txt_gz" : "data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz",
    "FastInsu_adjBMI_MAGIC_Europeans_txt_gz" : "data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz",
    "FastGlu_MAGIC_Europeans_txt_gz" : "data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz",
    "BMI_GIANT_2018_txt_gz" : "data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz",
    "MI_adjBMI_European_txt_gz" : "data/gwas/formatted/sumstats/hg38/MI_adjBMI_European/MI_adjBMI_European.txt.gz",
    "MAGIC_ISI_Model_2_AgeSexBMI_txt_txt_gz" : "data/gwas/formatted/sumstats/hg38/MAGIC_ISI_Model_2_AgeSexBMI.txt/MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz",
    "HDL_GLGC_Expanded_txt_gz" : "data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz",
    "TG_GLGC_Expanded_txt_gz": "data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz",
    "T2D_Mahajan_Europeans_txt_gz" : "data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz",
    "Type-2-Diabetes_Xue_2018_txt_gz" : "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz",
    "Type-2-Diabetes_Spracklen_2020_txt_gz" : "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz",
    "Type-2-Diabetes_Suzuki_2019_txt_gz" : "data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz"
}

# NOTE: If the ref SNP (top GWAS hit) is not actually in the eQTL file,
# then it currently won't work out for that locus

# First, load coloc info
coloc_tested_genes = {}
with open("output/post_coloc/de_genes/coloc_directionality.txt", "w") as w:
    w.write("hgnc\tgwas_file\teqtl_file\tchrom\tpos\tclpp_mod\tg_ea\tg_oa\tg_dir\te_ea\te_oa\te_dir\tsame_refs\tsame_effect_dirs\thigher_expression_higher_risk\n")
    with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split("\t")
            ensembl = data[11]
            hgnc = data[12]
            clpp_mod = float(data[14])
            eqtl_file = data[9]
            gwas_file = gwas_aliases[data[8].replace(".", "_")] 
            rsid = data[0]
            chrom = data[1]
            pos = data[2]

            # We have to just disregard sQTLs for now
            if "sQTL" in data[7]:
                continue

            # Is GWAS risk associated with increased or decreased expression?

            print gwas_file

            # Get GWAS SNP
            with gzip.open(gwas_file) as f:
                g_head = f.readline().strip().split()
            g_rsid_index = g_head.index("rsid")
            g_effect_allele_index = g_head.index("effect_allele")
            g_non_effect_allele_index = g_head.index("non_effect_allele")
            
            g_direction_index = g_head.index("effect_direction")
            gwas_data = subprocess.check_output("tabix {0} {1}:{2}-{2}".format(gwas_file, chrom, pos), shell=True).strip().split("\n")
            for gd in gwas_data:
                g_dat = gd.split("\t")
                if g_dat[g_rsid_index] != rsid:
                    continue
                else:
                    g_ea = g_dat[g_effect_allele_index]
                    g_oa = g_dat[g_non_effect_allele_index]
                    g_dir = g_dat[g_direction_index] == "+"
                    break

            # Get eQTL SNP
            with gzip.open(eqtl_file) as f:
                e_head = f.readline().strip().split()
            e_effect_allele_index = e_head.index("alt")
            e_non_effect_allele_index = e_head.index("ref")
            e_direction_index = e_head.index("beta")
            eqtl_data = subprocess.check_output("tabix {0} {1}:{2}-{2}".format(eqtl_file, chrom, pos), shell=True).strip().split("\n")
            for ed in eqtl_data:
                e_dat = ed.split("\t")
                if ensembl not in ed:
                    continue
                else:
                    e_ea = e_dat[e_effect_allele_index]
                    e_oa = e_dat[e_non_effect_allele_index]
                    e_dir = float(e_dat[e_direction_index]) > 0
                    break

            print g_ea, g_oa, g_dir, e_ea, e_oa, e_dir

            if not (((g_ea == e_ea) and (g_oa == e_oa)) or ((g_ea == e_oa) and (g_oa == e_ea))):
                print "discrepancy"
                continue 
            same_refs = ((g_ea == e_ea) and (g_oa == e_oa))
            print same_refs
            same_effect_dirs = (e_dir == g_dir)
            print same_effect_dirs
            higher_expression_higher_risk = (same_refs == same_effect_dirs)
            print higher_expression_higher_risk

            if ensembl not in coloc_tested_genes:
                coloc_tested_genes[ensembl] = []

            coloc_genes = [hgnc, gwas_file, eqtl_file, chrom, pos, clpp_mod, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir, same_refs, same_effect_dirs, same_refs == same_effect_dirs]


            w.write("{0}\n".format("\t".join([str(cg) for cg in coloc_genes])))
                
