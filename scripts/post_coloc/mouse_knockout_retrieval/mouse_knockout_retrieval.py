# TODO: For genes missed because weren't in the human to mouse mapping, maybe there's another
# way I can find orthologues for these. See what some of them are first though

# TODO: Filter to only ones in one of Ivan's categories of interest

# TODO: do it for IMPC too

# TODO: do it for Bruna's DE genes list

# TODO: enrichment test from binomial or hypergeometric test

# TODO: merge with coloc tables, integrate it into the pipeline


import requests
import sys
import copy

coloc_results_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt"

impc_geno_to_pheno_file = "data/mouse_knockouts/ALL_genotype_phenotype.csv.gz"
human_mouse_map_file = "data/mouse_knockouts/HMD_HumanPhenotype.rpt"
mp_ontology_file = "data/mouse_knockouts/VOC_MammalianPhenotype.rpt"
mp_pheno_tree_file = "data/mouse_knockouts/MPheno_OBO.ontology"
out_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/mouse_knockout_phenos.txt"

required_parents = ["MP:0010769", # mortality / aging -- contains lethal phenotype
                    "MP:0002873", # normal phenotype / no abnormality
                    "MP:0005375", # adipose
                    "MP:0005379", # endocrine / exocrine
                    "MP:0005376", # homeostasis
                    "MP:0005370", # liver / biliary
                    "MP:0005369"] # muscle

def main():

    mp_phenos = load_ontology(mp_ontology_file)
    mp_tree = load_pheno_tree(mp_pheno_tree_file)

    mp_phenos = filter_phenos(mp_phenos, mp_tree, required_parents)
    
    gene_set = get_genes_from_coloc(coloc_results_file)
    print "Number genes tested in coloc:"
    print len(list(gene_set))

    # Get map of human genes to mouse genes
    human_to_mouse = get_key_mouse_genes(gene_set, human_mouse_map_file)

    print "Number genes with detected mouse orthologue:"
    print len(human_to_mouse.keys())

    # Get phenotypes corresponding to the gene knockouts in MGI
    gene_phenos = get_gene_phenos(human_to_mouse)
    
    print "Number genes with at least one knockout mouse in MGI:"
    print len([gp for gp in gene_phenos if len(gene_phenos[gp]) != 0])

    # Write an output file mapping human genes to the corresponding mouse phenotypes
    write_output_file(gene_set, mp_phenos, gene_phenos, human_to_mouse, out_file)

def filter_phenos(mp_phenos, mp_tree, required_parents):
    filtered_phenos = copy.deepcopy(mp_phenos)

    for pheno in mp_phenos:
        valid = False
        current_pheno = pheno
        while current_pheno != "":
            if current_pheno in required_parents:
                valid = True
                break
            current_pheno = mp_tree.get(current_pheno, "")
        if not valid:
            del filtered_phenos[pheno]

    return filtered_phenos

def write_output_file(gene_set, mp_phenos, gene_phenos, human_to_mouse, out_file):
    with open(out_file, "w") as w:
        w.write("human_gene\tmouse_gene\tphenotypes\n")
        for gene in list(gene_set):
            if gene in gene_phenos:
                if len(gene_phenos[gene]) == 0:
                    pheno_formatted = "no_knockout_mice"
                else:
                    pheno_formatted = ",".join([pheno + ":" + mp_phenos[pheno].strip().replace(" ", "_") for pheno in gene_phenos[gene] if pheno in mp_phenos])
                if pheno_formatted == "":
                    pheno_formatted = "no_relevant_phenotypes_found"
            else:
                pheno_formatted = "no_mouse_ortholog_found"
            w.write("{0}\t{1}\t{2}\n".format(gene, human_to_mouse.get(gene, "NA"), pheno_formatted))

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

def get_key_mouse_genes(human_gene_set, human_mouse_map_file):
    human_to_mouse = {}
    with open(human_mouse_map_file) as f:
        for line in f:
            data = line.strip().split("\t")
            human = data[0]
            mouse = data[5]
            if human in human_gene_set:
                human_to_mouse[human] = mouse

    return human_to_mouse


if __name__ == "__main__":
    main()
