import glob
import subprocess
import gzip
import os

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

lipid_traits = [gwas_aliases["HDL_GLGC_Expanded_txt_gz"],
                gwas_aliases["TG_GLGC_Expanded_txt_gz"]]
whr_traits = [gwas_aliases["WHR-adj-BMI_GIANT_2018_txt_gz"]]
glucose_traits = [gwas_aliases["T2D_Mahajan_Europeans_txt_gz"],
                    gwas_aliases["Type-2-Diabetes_Xue_2018_txt_gz"],
                    gwas_aliases["Type-2-Diabetes_Spracklen_2020_txt_gz"],
                    gwas_aliases["Type-2-Diabetes_Suzuki_2019_txt_gz"],
                    gwas_aliases["MAGIC_ISI_Model_2_AgeSexBMI_txt_txt_gz"],
                    gwas_aliases["FastGlu_MAGIC_Europeans_txt_gz"],
                    gwas_aliases["FastInsu_adjBMI_MAGIC_Europeans_txt_gz"]]
                # MI not included because directionality not clear to me right now

if os.path.isfile("output/post_coloc/de_genes/not_found.txt"):
    os.remove("output/post_coloc/de_genes/not_found.txt")

# First, load coloc info

hgnc_to_ensembl_map = {}

coloced_genes = {}
liver_lipids_coloced_genes = {}
adipose_whr_coloced_genes = {}
muscle_glucose_coloced_genes = {}
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

        ## We have to just disregard sQTLs for now
        #if "sQTL" in data[7]:
        #    continue

        if clpp_mod < 0.35:
            continue

        # Is GWAS risk associated with increased or decreased expression?

        print gwas_file
        print eqtl_file
        # Get GWAS SNP
        with gzip.open(gwas_file) as f:
            g_head = f.readline().strip().split()
        g_rsid_index = g_head.index("rsid")
        g_effect_allele_index = g_head.index("effect_allele")
        g_non_effect_allele_index = g_head.index("non_effect_allele")
        if "effect_direction" in g_head:
            g_direction_index = g_head.index("effect_direction")
        gwas_data = subprocess.check_output("tabix {0} {1}:{2}-{2}".format(gwas_file, chrom, pos), shell=True).strip().split("\n")
        found = False
        for gd in gwas_data:
            g_dat = gd.split("\t")
            if g_dat[g_rsid_index] != rsid:
                continue
            else:
                g_ea = g_dat[g_effect_allele_index]
                g_oa = g_dat[g_non_effect_allele_index]
                g_dir = g_dat[g_direction_index] == "+"
                found = True
                break

        if not found:
            with open("output/post_coloc/de_genes/not_found.txt", "a") as a:
                a.write("{0}\t{1}\t{2}\n".format(hgnc, ensembl, gwas_file))
            continue

        # Get eQTL SNP
        with gzip.open(eqtl_file) as f:
            e_head = f.readline().strip().split()
        e_effect_allele_index = e_head.index("alt")
        e_non_effect_allele_index = e_head.index("ref")
        e_direction_index = e_head.index("beta")
        eqtl_data = subprocess.check_output("tabix {0} {1}:{2}-{2}".format(eqtl_file, chrom, pos), shell=True).strip().split("\n")
        found = False
        for ed in eqtl_data:
            e_dat = ed.split("\t")
            if ensembl not in ed:
                continue
            else:
                e_ea = e_dat[e_effect_allele_index]
                e_oa = e_dat[e_non_effect_allele_index]
                e_dir = float(e_dat[e_direction_index]) > 0
                found = True
                break

        if not found:
            with open("output/post_coloc/de_genes/not_found.txt", "a") as a:
                a.write("{0}\t{1}\t{2}\n".format(hgnc, ensembl, gwas_file))
            continue

        print rsid, chrom, pos, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir

        if not (((g_ea == e_ea) and (g_oa == e_oa)) or ((g_ea == e_oa) and (g_oa == e_ea))):
            print "discrepancy"
            continue 
        same_refs = ((g_ea == e_ea) and (g_oa == e_oa))
        print same_refs
        same_effect_dirs = (e_dir == g_dir)
        print same_effect_dirs
        higher_expression_higher_risk = (same_refs == same_effect_dirs)
        print higher_expression_higher_risk

        if ensembl not in coloced_genes:
            coloced_genes[ensembl] = []
            hgnc_to_ensembl_map[hgnc] = ensembl

        coloced_genes[ensembl].append((hgnc, gwas_file, eqtl_file, chrom, pos, clpp_mod, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir, same_refs, same_effect_dirs, higher_expression_higher_risk))

        if "Liver" in eqtl_file and gwas_file in lipid_traits:
            if ensembl not in liver_lipids_coloced_genes:
                liver_lipids_coloced_genes[ensembl] = []

            liver_lipids_coloced_genes[ensembl].append((hgnc, gwas_file, eqtl_file, chrom, pos, clpp_mod, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir, same_refs, same_effect_dirs, higher_expression_higher_risk))

        if "Adipose" in eqtl_file and gwas_file in whr_traits:
            if ensembl not in adipose_whr_coloced_genes:
                adipose_whr_coloced_genes[ensembl] = []

            adipose_whr_coloced_genes[ensembl].append((hgnc, gwas_file, eqtl_file, chrom, pos, clpp_mod, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir, same_refs, same_effect_dirs, higher_expression_higher_risk))

        if "Muscle" in eqtl_file and gwas_file in glucose_traits:
            if ensembl not in muscle_glucose_coloced_genes:
                muscle_glucose_coloced_genes[ensembl] = []

            muscle_glucose_coloced_genes[ensembl].append((hgnc, gwas_file, eqtl_file, chrom, pos, clpp_mod, g_ea, g_oa, g_dir, e_ea, e_oa, e_dir, same_refs, same_effect_dirs, higher_expression_higher_risk))


# Also load single-coloc genes
# NOTE: Currently this might miss a few -- I should really cross-reference with ENSGs to be sure
single_coloc_genes = set([])
with open("data/curated_gene_sets/single_coloc_genes_updated.txt") as f:
    for line in f:
        data = line.strip()
        if data in hgnc_to_ensembl_map:
            print data
            single_coloc_genes.add(hgnc_to_ensembl_map[data])

single_coloc_set = {scg: coloced_genes[scg] for scg in coloced_genes if scg in single_coloc_genes}

def annotate_gene_set(gene_set, gene_set_name):

    # For each perturbation-tissue combo...
    perturb_tissues = glob.glob("data/de_genes/*/*.txt")
    with open("output/post_coloc/de_genes/perturbation_by_{0}_with_sqtl.txt".format(gene_set_name), "w") as w:
        w.write("pert_tissue\tperturbation\tpert_direction\tgene\thgnc\tgwas_file\teqtl_file\tchrom\tpos\tclpp_mod\tg_ea\tg_oa\tg_dir\te_ea\te_oa\te_dir\tsame_refs\tsame_effect_dirs\thigher_expression_higher_risk\tperturbation_increases_risk\n")
        for pt in perturb_tissues:

            with open(pt) as f:

                f.readline()

                for line in f:
                    data = line.strip().split()
                   
                    if len(data) < 9:
                        continue

                    if data[9] != "TRUE":
                        continue

                    # For each that has actual coloc...
                    if data[2] in gene_set:
                       
                        direction = float(data[5]) > 0

                        for cg in gene_set[data[2]]:

                            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(data[0], data[1], direction, data[2], "\t".join([str(c) for c in cg]), direction == cg[-1]))


annotate_gene_set(coloced_genes, "coloc")
annotate_gene_set(single_coloc_set, "single_coloc")
annotate_gene_set(liver_lipids_coloced_genes, "liver_lipids_coloc")
annotate_gene_set(adipose_whr_coloced_genes, "adipose_whr_coloc")
annotate_gene_set(muscle_glucose_coloced_genes, "muscle_glucose_coloc")

