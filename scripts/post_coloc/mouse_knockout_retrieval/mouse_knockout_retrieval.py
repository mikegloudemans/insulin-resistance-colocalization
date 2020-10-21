### Lower priority

# TODO: For genes missed because weren't in the human to mouse mapping, maybe there's another
# way I can find orthologues for these. See what some of them are first though
#
# Status: I checked this and I'd estimate only about 10% of cases with no found ortholog
# are ones where I actually missed a real one (e.g. AARS is properly listed as AARS1 in the lookup table).
# Solution would be to map from ENSG directly to mouse gene, if I can find a table for that,
# since ENSG aren't so ambiguous. Maybe not top priority right now though.

# TODO: merge with coloc and DE tables, integrate it into the fuller pipeline

### Top priority

# TODO: do it for Bruna's DE genes list

# TODO: enrichment test from binomial or hypergeometric test


import requests
import sys
import copy
import gzip

coloc_results_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt"
de_results_file = "data/mouse_knockouts/IR_GxE_DE_and_Matched_nonDE_genes.tsv"

impc_gene_to_pheno_file = "data/mouse_knockouts/ALL_genotype_phenotype.csv.gz"
human_mouse_map_file = "data/mouse_knockouts/HMD_HumanPhenotype.rpt"
mp_ontology_file = "data/mouse_knockouts/VOC_MammalianPhenotype.rpt"
mp_pheno_tree_file = "data/mouse_knockouts/MPheno_OBO.ontology"
out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockouts/{0}_mouse_knockout_phenos_mgi.txt"
impc_out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockouts/{0}_mouse_knockout_phenos_impc.txt"

required_parents = ["MP:0010769", # mortality / aging -- contains lethal phenotype
                    "MP:0002873", # normal phenotype / no abnormality
                    "MP:0005375", # adipose
                    "MP:0005379", # endocrine / exocrine
                    "MP:0005376", # homeostasis
                    "MP:0005370", # liver / biliary
                    "MP:0005369"] # muscle

def main():
  
    coloc_gene_set = get_genes_from_coloc(coloc_results_file)
    de_gene_set = get_genes_from_de(de_results_file)

    retrieve_and_write_knockouts(coloc_gene_set, "coloc")    
    #retrieve_and_write_knockouts(de_gene_set, "de")    

def retrieve_and_write_knockouts(gene_set, prefix):

    mp_phenos = load_ontology(mp_ontology_file)
    mp_tree = load_pheno_tree(mp_pheno_tree_file)

    mp_phenos, mp_parents = filter_phenos(mp_phenos, mp_tree, required_parents)
  
    print "-----------------------"

    print "Set:", prefix

    # Get map of human genes to mouse genes
    human_to_mouse = get_key_mouse_genes(gene_set, human_mouse_map_file)
    
    print "Number genes tested:"
    print len(list(gene_set))

    print "Number genes with detected mouse orthologue:"
    print len(human_to_mouse.keys())

    # Get phenotypes corresponding to the gene knockouts in MGI
    #impc_gene_phenos = get_impc_gene_phenos(human_to_mouse, impc_gene_to_pheno_file)
    gene_phenos = get_gene_phenos(human_to_mouse)
    
    print "Number genes with at least one knockout mouse in MGI:"
    print len([gp for gp in gene_phenos if len(gene_phenos[gp]) != 0])

    print "Number genes with at least one knockout mouse in IMPC:"
    #print len([gp for gp in impc_gene_phenos if len(impc_gene_phenos[gp]) != 0])

    # Write an output file mapping human genes to the corresponding mouse phenotypes
    write_output_file(gene_set, mp_phenos, mp_parents, gene_phenos, human_to_mouse, out_file.format(prefix))
    #write_output_file(gene_set, mp_phenos, mp_parents, impc_gene_phenos, human_to_mouse, impc_out_file.format(prefix))

def filter_phenos(mp_phenos, mp_tree, required_parents):
    filtered_phenos = copy.deepcopy(mp_phenos)

    parents = {}

    for pheno in mp_phenos:
        valid = False
        current_pheno = pheno
        while current_pheno != "":
            if current_pheno in required_parents:
                parents[pheno] = current_pheno
                valid = True
                break
            current_pheno = mp_tree.get(current_pheno, "")
        if not valid:
            del filtered_phenos[pheno]

    return filtered_phenos, parents

def write_output_file(gene_set, mp_phenos, mp_parents, gene_phenos, human_to_mouse, out_file):
    with open(out_file, "w") as w:
        w.write("human_gene\tmouse_gene\tphenotypes\tparents\n")
        for gene in sorted(list(gene_set)):
            if gene in gene_phenos:
                if len(gene_phenos[gene]) == 0:
                    pheno_formatted = "no_knockout_mice"
                    parents = "NA"
                else:
                    pheno_formatted = ",".join([pheno + ":" + mp_phenos[pheno].strip().replace(" ", "_").replace(",", "") for pheno in gene_phenos[gene] if pheno in mp_phenos])
                    parents = list(set([mp_parents[pheno] + ":" + mp_phenos[mp_parents[pheno]].strip().replace(" ", "_").replace(",", "") for pheno in gene_phenos[gene] if pheno in mp_phenos]))
                if pheno_formatted == "":
                    pheno_formatted = "no_relevant_phenotypes_found"
            else:
                pheno_formatted = "no_mouse_ortholog_found"
            w.write("{0}\t{1}\t{2}\t{3}\n".format(gene, human_to_mouse.get(gene, "NA"), pheno_formatted, parents))

def load_pheno_tree(mp_pheno_tree_file):
    pheno_tree = {}
    with open(mp_pheno_tree_file) as f:
        current_term = ""
        for line in f:
            if line.startswith("id:"):
                current_term = line.strip().replace("id: ", "")
            if line.startswith("is_a:"):
                pheno_tree[current_term] = line.strip().replace("is_a: ", "").split("!")[0].strip()

    return pheno_tree


def load_ontology(mp_ontology_file):
    mp_phenos = {}
    with open(mp_ontology_file) as f:
        for line in f:
            data = line.strip().split("\t")
            mp_phenos[data[0]] = data[1]

    return mp_phenos


def get_gene_phenos(human_to_mouse):

    gene_phenos = {}

    for hgene in human_to_mouse.keys():
        mgene = human_to_mouse[hgene]

        r = requests.get("http://www.informatics.jax.org/marker/phenotypes/{0}".format(mgene))
        text = r.text

        terms = text.split("http://www.informatics.jax.org/vocab/mp_ontology/")[1:]
        terms = [t.split('\">')[0] for t in terms]
        gene_phenos[hgene] = terms

    return gene_phenos

def get_impc_gene_phenos(human_to_mouse, impc_gene_to_pheno_file):
    
    gene_phenos = {}
    # Reverse the dictionary for this step
    mouse_to_human = {v: k for k, v in human_to_mouse.iteritems()}

    with gzip.open(impc_gene_to_pheno_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split(",")
            mgene = data[0]
            pheno = [d for d in data if "MP:" in d]

            if mgene not in mouse_to_human:
                continue
            hgene = mouse_to_human[mgene]

            if hgene not in gene_phenos:
                gene_phenos[hgene] = []

            if pheno not in gene_phenos[hgene]:
                gene_phenos[hgene].extend(pheno)

    return gene_phenos


def get_genes_from_coloc(coloc_results_file):
    
    gene_set = set([])

    # First, load the coloc results
    # and get all genes from them
    with open(coloc_results_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split("\t")
            gene = data[12]
            gene_set.add(gene)

    return gene_set

def get_genes_from_de(de_results_file):
    
    # For now, it doesn't matter if they're actually DE
    # or if they're null controls, because we're just fetching
    # the annotations.

    de_gene_set = set([])

    with open(de_results_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            gene = data[2]
            de_gene_set.add(gene)

    return de_gene_set

def get_key_mouse_genes(human_gene_set, human_mouse_map_file):
    human_to_mouse = {}
    with open(human_mouse_map_file) as f:
        for line in f:
            data = line.strip().split("\t")
            human = data[0]
            mouse = data[5].strip()
            if human in human_gene_set:
                human_to_mouse[human] = mouse

    return human_to_mouse


if __name__ == "__main__":
    main()
