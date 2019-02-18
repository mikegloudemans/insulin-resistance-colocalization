import gzip

# Load rsid mapping
rsids = {}
with gzip.open("/users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz") as f:
    for line in f:
        data = line.strip().split()
        rsids[(data[0].replace("chr", ""), int(data[1]))] = data[2]

# List all studies to browse through
# If possible, use recently formatted ones since they fit the coloc pipeline best?
studies = [
            "/users/mgloud/projects/gwas/data/prepared/FastInsu_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/FastGlu_MAGIC_Europeans_AllSNPs.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/TG_GCLC_Mixed.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/LDL_GCLC_Mixed.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/HDL_GCLC_Mixed.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/GENESIS/formatted/GWAS_GENESIS_MI.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/T2D_DIAGRAM_European.prepared.txt.gz",
            "/users/mgloud/projects/gwas/data/prepared/BMI_GIANT_Europeans_AllSNPs.prepared.txt.gz"
        ]

hits = []
# For each study
for study in studies:
    print study
    with gzip.open(study) as f:
        header = f.readline().strip().split()
        chr_index = header.index("chr")
        pos_index = header.index("snp_pos")
        pval_index = header.index("pvalue")
        #rsid_index = header.index("rsid")
        for line in f:
            data = line.strip().split()
            try:
                pval = float(data[pval_index])
            except:
                print line
                continue
            if study == "/users/mgloud/projects/gwas/data/GENESIS/formatted/GWAS_GENESIS_MI.txt.gz":
                if pval > 2e-6:
                    continue
            elif pval > 5e-8:
                continue
            if "hap" in data[chr_index]:
                continue
            hits.append((study, data[chr_index].replace('chr', ''), int(data[pos_index]), pval))

print hits[:10]

# Make list of all SNPs significant in at least one study
hit_map = {}
for hit in hits:
    if (hit[1], hit[2]) not in hit_map:
        hit_map[(hit[1], hit[2])] = {"studies": [], "gwas_pvals": [], "genes": [], "eqtl_pvals": []}

    hit_map[(hit[1], hit[2])]["studies"].append(hit[0].split("/")[-1])
    hit_map[(hit[1], hit[2])]["gwas_pvals"].append(hit[3])


# Get eQTL p-values for every gene, for each SNP

with gzip.open("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/all_associations/Liver.allpairs.txt.gz") as f:
    f.readline()
    for line in f:
        data = line.strip().split()
        snpdata = data[1].split("_")
        if (snpdata[0], int(snpdata[1])) in hit_map:
            hit_map[(snpdata[0], int(snpdata[1]))]["genes"].append(data[0])
            hit_map[(snpdata[0], int(snpdata[1]))]["eqtl_pvals"].append(data[6])

with open("/users/mgloud/projects/gwas/output/insulin_snps_to_test.txt", "w") as w:
    w.write("rsid\tchrom\tposition\tgwas_traits\tgwas_pvals\teqtl_genes\teqtl_pvals\n")
    for k in sorted(hit_map.keys()):
        if (k[0], k[1]) in rsids:
            rsid = rsids[(k[0], k[1])]
        else:
            rsid = "rsid_unknown"
        w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(rsid, k[0], k[1], ",".join(hit_map[k]["studies"]), ",".join([str(g) for g in hit_map[k]["gwas_pvals"]]), ",".join(hit_map[k]["genes"]), ",".join(hit_map[k]["eqtl_pvals"])))

