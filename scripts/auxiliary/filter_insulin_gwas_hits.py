# 08/2019 -- I can't remember what the point of this was before, but I don't really think it's
# useful anymore

import operator

# Load aill GWAS p-values
gwas = {}
every_snp = []
with open("/users/mgloud/projects/gwas/output/insulin_snps_to_test.txt") as f:
    f.readline()
    for line in f:
        data = line.strip().split()
        traits = data[3].strip().split(",")
        pvals = [float(p) for p in data[4].strip().split(",")]
        
        assert len(traits) == len(pvals)

        for i in range(len(traits)):
            if traits[i] not in gwas:
                gwas[traits[i]] = []

            gwas[traits[i]].append((data[0], pvals[i], int(data[1]), int(data[2])))
            every_snp.append((data[0], pvals[i], int(data[1]), int(data[2])))

# For each trait, get the set of independent loci for which there is
# a genome-wide significant hit.
all_hits = []
for trait in gwas:
    top_hits = []
    ordered = sorted(gwas[trait], key=operator.itemgetter(1))

    for o in ordered:
        good = True
        for th in top_hits:
            if o[2] == th[2] and abs(o[3] - th[3]) < 1000000:
                good = False
                break
        if good:
            top_hits.append(o)

    all_hits.extend(top_hits)
all_rsids = set([ah[0] for ah in all_hits if ah[0] != "rsid_unknown"])

# Write a new file containing only the top hits for each trait
with open("/users/mgloud/projects/insulin_resistance/output/filtered_insulin_snps_to_test.txt", "w") as w:
    with open("/users/mgloud/projects/gwas/output/insulin_snps_to_test.txt") as f:
        w.write(f.readline())
        for line in f:
            data = line.strip().split()
            if data[0] in all_rsids:
                w.write(line)

# Now just the key info to run coloc
with open("/users/mgloud/projects/insulin_resistance/output/coloc_ready_insulin_snp_list.txt", "w") as w:
    with open("/users/mgloud/projects/gwas/output/insulin_snps_to_test.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            if data[0] in all_rsids:
                w.write("\t".join(data[1:3]) + "\n")

# Get the set of independent loci for which there is a genome-wide
# significant hit for ANY trait
top_hits = []
ordered = sorted(every_snp, key=operator.itemgetter(1))
for o in ordered:
    good = True
    for th in top_hits:
        if o[2] == th[2] and abs(o[3] - th[3]) < 1000000:
            good = False
            break
    if good:
        top_hits.append(o)
all_rsids = set([ah[0] for ah in top_hits if ah[0] != "rsid_unknown"])

with open("/users/mgloud/projects/insulin_resistance/output/filtered_clustered_insulin_snps_to_test.txt", "w") as w:
    with open("/users/mgloud/projects/gwas/output/insulin_snps_to_test.txt") as f:
        w.write(f.readline())
        for line in f:
            data = line.strip().split()
            if data[0] in all_rsids:
                w.write(line)
