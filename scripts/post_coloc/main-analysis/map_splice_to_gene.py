# Author: Mike Gloudemans
# Date: 5/25/2020
#
# Make a simple map from all hg38 splice sites to the most likely
# gene(s) containing them.
#
# Some splice sites may be assigned to more than one gene.
#

import gzip
import glob
import subprocess

def main():

    # Map each splice site to the associated gene or gene(s)
    splice_map = map_splice_junctions_to_genes()

    print splice_map.items()[:10]

    # Output a table showing the most probable mappings
    write_splice_map(splice_map, "output/splice_map/splice_map.txt")

def map_splice_junctions_to_genes():

    splice_map = {}

    tissues = glob.glob("data/sqtls/gtex_v8/*.gz")

    for tissue in tissues[:1]:
        print "Getting splice sites from:", tissue
        with gzip.open(tissue) as f:
            i = 0
            for line in f:
                i += 1
                if i > 100000:
                    break
                data = line.strip().split()
                junction, gene = data[0], data[1]
                if junction not in splice_map:
                    splice_map[junction] = set([])
                splice_map[junction].add(gene)

    return splice_map

def write_splice_map(splice_map, out_file):
    with open(out_file, "w") as w:
        for junction in splice_map:
            pass

if __name__ == "__main__":
    main()
