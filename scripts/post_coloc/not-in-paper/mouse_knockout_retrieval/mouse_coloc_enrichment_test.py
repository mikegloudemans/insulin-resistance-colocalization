from scipy import stats
import numpy as np

pli_file = "data/ExAC/release0.3/functional_gene_constraint/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt"

coloc_mgi_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_mouse_knockout_phenos_mgi.txt"
coloc_impc_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_mouse_knockout_phenos_impc.txt"
coloc_genes_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt"

coloc_mgi_out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_mgi_enrichments.txt"
coloc_impc_out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/coloc_impc_enrichments.txt"

# If we have fewer than this many genes annotated with a term, 
# don't even bother testing it
min_gene_limit = 0

def main():

    print 
    print "Testing MGI enrichment"
    test_enrichment(coloc_mgi_file, coloc_mgi_out_file, coloc_genes_file)
    print 
    print "Testing IMPC enrichment"
    test_enrichment(coloc_impc_file, coloc_impc_out_file, coloc_genes_file)

def load_pli():

    pli = {}

    with open(pli_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            gene = data[1]
            score = float(data[19])
            pli[gene] = score

    return pli

def test_enrichment(knockout_file, out_file, coloc_genes_file):

    # Get phenotypes for each gene
    gene_phenos = {}
    with open(knockout_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            gene = data[0]

            phenos = data[2].strip().split(",")

            gene_phenos[gene] = phenos

    # Make list of which genes are coloc genes and which are tested
    # but not coloc
    coloc_status = {}
    with open(coloc_genes_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            gene = data[12]
            clpp_mod = float(data[14])
            if (clpp_mod > 0.35):
                coloc_status[gene] = 1
            elif gene not in coloc_status:
                coloc_status[gene] = 0

    count_matrix = {}
    count_matrix["total_coloc_with_knockouts"] = 0
    count_matrix["total_control_with_knockouts"] = 0
    count_matrix["phenos"] = {}
        
    for gene in coloc_status:

        status = bool(coloc_status[gene])

        # Every gene should have SOME annotation in the phenotype list,
        # even if it's just "no knockout" or "no mouse ortholog"
        assert gene in gene_phenos

        # Don't worry about genes if there wasn't a mouse knockout for them
        if gene_phenos[gene][0] == "no_knockout_mice" or \
                gene_phenos[gene][0] == "no_mouse_ortholog_found":
            continue

        if status:
            count_matrix["total_coloc_with_knockouts"] += 1
        else:
            count_matrix["total_control_with_knockouts"] += 1

        phenos_hit = set([])

        for pheno in gene_phenos[gene]:
            # A pheno should only be counted once per gene
            if pheno in phenos_hit:
                continue
            phenos_hit.add(pheno)

            if pheno not in count_matrix["phenos"]: 
                count_matrix["phenos"][pheno] = {}
                count_matrix["phenos"][pheno]["coloc"] = 0
                count_matrix["phenos"][pheno]["control"] = 0

            if status:
                count_matrix["phenos"][pheno]["coloc"] += 1
            else:
                count_matrix["phenos"][pheno]["control"] += 1
                
    # Write enrichment stats to an output file:
    with open(out_file, "w") as w:

        w.write("mouse_phenotype\tnum_coloc_genes\tnum_coloc_ko_genes_without_term\tfraction_coloc_genes\tnum_control_genes\tnum_control_ko_genes_without_term\tfraction_control_genes\tfishers_test_pvalue\n")
        
        for pheno in count_matrix["phenos"]:
            if count_matrix["phenos"][pheno]["coloc"] + \
                    count_matrix["phenos"][pheno]["control"] < min_gene_limit:
                continue
            cont_table = [[count_matrix["phenos"][pheno]["coloc"], \
                        count_matrix["total_coloc_with_knockouts"] - count_matrix["phenos"][pheno]["coloc"]  ], \
                    [count_matrix["phenos"][pheno]["control"],\
                        count_matrix["total_control_with_knockouts"] - count_matrix["phenos"][pheno]["control"]]]
            pvalue = stats.fisher_exact(cont_table)[1]
            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(\
                    pheno, \
                            count_matrix["phenos"][pheno]["coloc"], \
                            count_matrix["total_coloc_with_knockouts"] - count_matrix["phenos"][pheno]["coloc"], \
                            count_matrix["phenos"][pheno]["coloc"] * 1.0 / count_matrix["total_coloc_with_knockouts"],\
                            count_matrix["phenos"][pheno]["control"],\
                            count_matrix["total_control_with_knockouts"] - count_matrix["phenos"][pheno]["control"], \
                            count_matrix["phenos"][pheno]["control"] * 1.0 / count_matrix["total_control_with_knockouts"],\
                            pvalue))

    # Compare PLI scores while we're at it

    pli = load_pli()

    coloc_plis = []
    non_coloc_plis = []
    for gene in coloc_status:
        
        if gene not in pli:
            continue

        status = bool(coloc_status[gene])

        if status:
            coloc_plis.append(pli[gene])
        else:
            non_coloc_plis.append(pli[gene])

    print "coloc pLI distribution:"
    print [round(n*100,2) for n in np.quantile(coloc_plis, [0.01 * x for x in range(0,101,5)])]
    print np.median(coloc_plis)
    print
    print "non-coloc pLI distribution:"
    print [round(n*100,2) for n in np.quantile(non_coloc_plis, [0.01 * x for x in range(0,101,5)])]
    print np.median(non_coloc_plis)


if __name__ == "__main__":
    main()
