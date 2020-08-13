from scipy import stats

de_genes_file = "data/mouse_knockouts/IR_GxE_DE_and_Matched_nonDE_genes.tsv"
de_mgi_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/de_mouse_knockout_phenos_mgi.txt"
mgi_out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/de_mgi_enrichments.txt"

de_impc_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/de_mouse_knockout_phenos_impc.txt"
impc_out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/de_impc_enrichments.txt"

# If we have fewer than this many genes annotated with a term, 
# don't even bother testing it
min_gene_limit = 0

def main():

    test_enrichment(de_mgi_file, mgi_out_file, de_genes_file)
    test_enrichment(de_impc_file, impc_out_file, de_genes_file)

def test_enrichment(knockout_file, out_file, de_genes_file):

    # Get phenotypes for each gene
    gene_phenos = {}
    with open(knockout_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            gene = data[0]

            phenos = data[2].strip().split(",")

            gene_phenos[gene] = phenos

    # Count up, for each cell type and condition,
    # the number of DE and non-DE genes for which
    # each phenotype was detected in mice

    count_matrix = {}
    with open(de_genes_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()

            cell = data[0]
            condition = data[1]
            gene = data[2]
            de_status = bool(int(data[3]))

            if cell not in count_matrix:
                count_matrix[cell] = {}
            if condition not in count_matrix[cell]:
                count_matrix[cell][condition] = {}
                count_matrix[cell][condition]["total_de_with_knockouts"] = 0
                count_matrix[cell][condition]["total_control_with_knockouts"] = 0
                count_matrix[cell][condition]["phenos"] = {}

            # Every gene should have SOME annotation in the phenotype list,
            # even if it's just "no knockout" or "no mouse ortholog"
            assert gene in gene_phenos

            # Don't worry about genes if there wasn't a mouse knockout for them
            if gene_phenos[gene][0] == "no_knockout_mice" or \
                    gene_phenos[gene][0] == "no_mouse_ortholog_found":
                continue

            if de_status:
                count_matrix[cell][condition]["total_de_with_knockouts"] += 1
            else:
                count_matrix[cell][condition]["total_control_with_knockouts"] += 1

            phenos_hit = set([])

            for pheno in gene_phenos[gene]:
                # A pheno should only be counted once per gene
                if pheno in phenos_hit:
                    continue
                phenos_hit.add(pheno)

                if pheno not in count_matrix[cell][condition]["phenos"]: 
                    count_matrix[cell][condition]["phenos"][pheno] = {}
                    count_matrix[cell][condition]["phenos"][pheno]["de"] = 0
                    count_matrix[cell][condition]["phenos"][pheno]["control"] = 0

                if de_status:
                    count_matrix[cell][condition]["phenos"][pheno]["de"] += 1
                else:
                    count_matrix[cell][condition]["phenos"][pheno]["control"] += 1
                    
    # Write enrichment stats to an output file:
    with open(out_file, "w") as w:

        w.write("cell\tcondition\tmouse_phenotype\tnum_de_genes\tnum_de_ko_genes_without_term\tfraction_de_genes\tnum_control_genes\tnum_control_ko_genes_without_term\tfraction_control_genes\tfishers_test_pvalue\n")
        
        for cell in count_matrix:
            for condition in count_matrix[cell]:
                for pheno in count_matrix[cell][condition]["phenos"]:
                    if count_matrix[cell][condition]["phenos"][pheno]["de"] + \
                            count_matrix[cell][condition]["phenos"][pheno]["control"] < min_gene_limit:
                        continue
                    cont_table = [[count_matrix[cell][condition]["phenos"][pheno]["de"], \
                                count_matrix[cell][condition]["total_de_with_knockouts"] - count_matrix[cell][condition]["phenos"][pheno]["de"]  ], \
                            [count_matrix[cell][condition]["phenos"][pheno]["control"],\
                                count_matrix[cell][condition]["total_control_with_knockouts"] - count_matrix[cell][condition]["phenos"][pheno]["control"]]]
                    pvalue =  stats.fisher_exact(cont_table)[1]
                    w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(\
                            cell, condition, pheno, \
                                    count_matrix[cell][condition]["phenos"][pheno]["de"], \
                                    count_matrix[cell][condition]["total_de_with_knockouts"] - count_matrix[cell][condition]["phenos"][pheno]["de"], \
                                    count_matrix[cell][condition]["phenos"][pheno]["de"] * 1.0 / count_matrix[cell][condition]["total_de_with_knockouts"],\
                                    count_matrix[cell][condition]["phenos"][pheno]["control"],\
                                    count_matrix[cell][condition]["total_control_with_knockouts"] - count_matrix[cell][condition]["phenos"][pheno]["control"], \
                                    count_matrix[cell][condition]["phenos"][pheno]["control"] * 1.0 / count_matrix[cell][condition]["total_control_with_knockouts"],\
                                    pvalue))

if __name__ == "__main__":
    main()
