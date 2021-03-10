sqtl_ensembl = set([])
eqtl_ensembl = set([])
sqtl_locus = set([])
eqtl_locus = set([])
sqtl_unique_ensembl = set([])
eqtl_unique_ensembl = set([])
sqtl_unique_locus = set([])
eqtl_unique_locus = set([])

with open ("output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/clpp_results_categorized_2020-05-11.txt") as f:
    f.readline()
    for line in f:
        data = line.strip().split()
        clpp_mod = float(data[14])

        if clpp_mod < 0.35:
            continue

        if "sQTL" in data[7]:
            sqtl_ensembl.add(data[11])
            sqtl_locus.add(data[3])
            if data[19] in ["loci_2", "loci_0.1"]:
                sqtl_unique_ensembl.add(data[11])
                sqtl_unique_locus.add(data[3])
        elif "eQTL" in data[7]:
            eqtl_ensembl.add(data[11])
            eqtl_locus.add(data[3])
            if data[19] in ["loci_2", "loci_0.1"]:
                eqtl_unique_ensembl.add(data[11])
                eqtl_unique_locus.add(data[3])
        else:
            print "something's wrong"

print "loci"
print(len(list(sqtl_locus)))
print(len(list(eqtl_locus)))
print(len(list(eqtl_locus.intersection(sqtl_locus))))

print "unique loci"
print(len(list(sqtl_unique_locus)))
print(len(list(eqtl_unique_locus)))
print(len(list(eqtl_unique_locus.intersection(sqtl_unique_locus))))

print "ensembl"
print(len(list(sqtl_ensembl)))
print(len(list(eqtl_ensembl)))
print(len(list(eqtl_ensembl.intersection(sqtl_ensembl))))

print "unique ensembl"
print(len(list(sqtl_unique_ensembl)))
print(len(list(eqtl_unique_ensembl)))
print(len(list(eqtl_unique_ensembl.intersection(sqtl_unique_ensembl))))
