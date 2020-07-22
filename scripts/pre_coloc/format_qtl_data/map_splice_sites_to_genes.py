import gzip
import glob

mappings = {}
tissues = glob.glob("data/sqtls/gtex_v8/*.sQTLs.txt.gz")
with gzip.open("data/sqtls/splice_to_gene_map.txt.gz", "w") as w:
    for tissue in tissues:
        print tissue
        w.write("splice_junction\tgene\n")
        with gzip.open(tissue) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                splice, gene = data[0], data[1].split(":")[1]
                if splice not in mappings:
                    mappings[splice] = []
                if gene not in mappings[splice]:
                    mappings[splice].append(gene)
                    w.write("{0}\t{1}\n".format(splice, gene))
