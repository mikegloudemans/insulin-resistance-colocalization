import gzip

# Get the lead SNP for the top coloc for each gene
lead_snp_map = {}
with open("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt") as f:
    f.readline()
    for line in f:
        data = line.strip().split("\t")
        chrom, pos = data[1], int(data[2])
        gene = data[12]
        ensg = data[11]
        clpp_mod = float(data[14])

        if gene not in lead_snp_map:
            lead_snp_map[gene] = ((chrom, pos), clpp_mod, ensg)
        elif lead_snp_map[gene][1] < clpp_mod:
            lead_snp_map[gene] = ((chrom, pos), clpp_mod, ensg)

# Load gencode exons
all_exons = []
with gzip.open("data/gencode/gencode.v31.annotation.gtf.gz") as f:
    for line in f:
        if line.startswith("#"):
            continue
        data = line.strip().split("\t")
        if data[2] != "exon":
            continue
        gene_id = data[8].split(".")[0].replace("gene_id \"", "").replace("\"", "")
        gene_name = data[8].split(" gene_name \"")[1].split("\";")[0]
        chrom = data[0]
        midpoint = (int(data[3]) + int(data[4])) / 2
        all_exons.append((chrom, midpoint, gene_id, gene_name))

# Get background set
gene_set = []
with open("data/curated_gene_sets/single_coloc_genes_updated.txt") as f:
    for line in f:
        gene = line.strip()

        # Using coloc results, figure out what the "seed" GWAS SNP was here for the top SNP
        if gene in lead_snp_map:
            print gene, lead_snp_map[gene][1]
            lead_snp = lead_snp_map[gene][0]
            lead_chrom = "chr" + lead_snp[0]
            lead_ensg = lead_snp_map[gene][2]
            lead_pos = lead_snp[1]
        else:
            print gene, "not found"
            continue

        # Then using GENCODE, get the closest gene for that one
        closest = (-1, 1000000000, -1)
        for ae in all_exons:
            if lead_chrom != ae[0]:
                continue
            if abs(lead_pos - ae[1]) < abs(lead_pos - closest[1]):
                closest = ae
        gene_set.append((gene, closest, lead_snp, lead_ensg))

collisions = 0

# Output background genes to a list then, for enrichment testing
with open("data/curated_gene_sets/single_coloc_genes_updated_matched_background.txt", "w") as w:
    w.write("coloc_gene\tcoloc_ensembl\tnearest_ensembl\tnearest_gene\tsnp\n")
    for gs in gene_set:
        print gs
        # Identical ENSG
        if gs[1][2].split(".")[0] == gs[0]:
            collisions += 1
        elif gs[0] == gs[1][3]:
            collisions += 1
        elif gs[0] == gs[1][3].replace("-AS1", "") or gs[0].replace("-AS1", "") == gs[1][3]:
            collisions += 1

        w.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gs[0], gs[3], gs[1][2], gs[1][3], "_".join([str(g) for g in gs[2]])))

print collisions
