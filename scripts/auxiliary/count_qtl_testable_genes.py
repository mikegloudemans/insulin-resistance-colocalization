import gzip
testable = set([])

eqtl_files = ["data/eqtls/gtex_v8/Pancreas.allpairs.txt.gz.eQTLs.txt.gz",
                "data/eqtls/gtex_v8/Liver.allpairs.txt.gz.eQTLs.txt.gz",
                "data/eqtls/gtex_v8/Muscle_Skeletal.allpairs.txt.gz.eQTLs.txt.gz",
                "data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz",
                "data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz"]
sqtl_files = ["data/sqtls/gtex_v8/Pancreas.sQTLs.txt.gz",
                "data/sqtls/gtex_v8/Muscle_Skeletal.sQTLs.txt.gz",
                "data/sqtls/gtex_v8/Liver.sQTLs.txt.gz",
                "data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz",
                "data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz"]

for ef in eqtl_files:
    print ef
    with gzip.open(ef) as f:
        f.readline()
        for line in f:
            gene = line.strip().split()[0].split(".")[0]
            testable.add(gene)

for sf in sqtl_files:
    print sf
    with gzip.open(sf) as f:
        f.readline()
        for line in f:
            gene = line.strip().split()[1].split(":")[1].split(".")[0]
            testable.add(gene)

print len(list(testable)), "tested genes"
