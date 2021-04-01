coloc_set = set([])

with open("data/curated_gene_sets/single_coloc_genes_updated.txt") as f:
        for line in f:
            if line != "":
                coloc_set.add(line.strip())

with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockouts/single_coloc_mouse_summary.txt", "w") as w:
    w.write("Gene\tendocrine / exocrine\tadipose\thomeostasis\tmuscle\tliver / biliary\tsurvival\tnormal KO\tKO found, phenos not relevant\tno KO\tortholog not found\n")
    with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockouts/coloc_mouse_knockout_phenos_mgi.txt") as f:
        for line in f:
            if line.strip() == "":
                continue
            data = line.strip().split("\t")
            if data[0] in coloc_set:
                
                no_ko = "X" if data[2] == "no_knockout_mice" else ""
                no_mouse_ortholog = "X" if data[2] == "no_mouse_ortholog_found" else ""
                no_pheno = "X" if data[2] == "no_relevant_phenotypes_found" else ""
                normal = "X" if "MP:0002873:normal_phenotype" in data[3] else ""
                endocrine = "X" if "MP:0005379:endocrine/exocrine_gland_phenotype" in data[3] else ""
                homeostasis = "X" if "MP:0005376:homeostasis/metabolism_phenotype" in data[3] else ""
                survival = "X" if "MP:0010769:abnormal_survival" in data[3] else ""
                muscle = "X" if "MP:0005369:muscle_phenotype" in data[3] else ""
                adipose = "X" if "MP:0005375:adipose_tissue_phenotype" in data[3] else ""
                biliary = "X" if "MP:0005370:liver/biliary_system_phenotype" in data[3] else ""

                w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(data[0],endocrine,adipose,homeostasis,muscle,biliary,survival,normal,no_pheno,no_ko,no_mouse_ortholog))

                coloc_set.remove(data[0])
                
    # Now write all the ones that we missed

    for cs in list(coloc_set):
        w.write("{0}\t\t\t\t\t\t\t\t\t\t\n".format(cs))
