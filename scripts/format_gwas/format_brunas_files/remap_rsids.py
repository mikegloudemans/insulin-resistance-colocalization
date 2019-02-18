# Map hg19 GWAS genomic coordinates to rsIDs, for use in Bosh's web application.
#
# rsIDs determined from 1Kgenomes.
#
# NOTE Warning: some rsID matchings may be incorrect, since we're not doing any testing to make sure the ref/alt allele are matched in GWAS and in 1kgenomes.
#

import sys
import subprocess
import gzip

file = sys.argv[1]
out = sys.argv[2]

trait_mode = False
with gzip.open(file) as f:
    header = f.readline().strip().split()

    snp_col = header.index("snp_pos") + 1
    chr_col = header.index("chr") + 1
    pvalue_col = header.index("pvalue") + 1
    
    if "trait" in header:
        trait_mode = True
        trait_col = header.index("trait") + 1


if not trait_mode:
    # Write key data to temporary file
    subprocess.check_call('''zcat {3} | tail -n +2 | awk '{{print ${0} "_" ${1} "\t" ${2} }}' | sed s/chr//g | sort -k1,1 > {4}_tmp'''.format(chr_col, snp_col, pvalue_col, file, out), shell=True)

    # Write header to final file
    subprocess.check_call('echo "chr\tpos\trsid\tpval" > {0}'.format(out), shell=True)

    # Write rsids to final file
    subprocess.check_call('''join /users/mgloud/projects/gwas/data/rsid_mapping_sorted.txt {0}_tmp | grep rs | grep -v \; | grep -v , | sed s/\ /\\\t/g | sed s/_/\\\t/g | sort -k1,1n -k2,2n >> {0}'''.format(out), shell=True)

    subprocess.check_call("rm {0}_tmp".format(out), shell=True)

else:
    # Write key data to temporary file
    subprocess.check_call('''zcat {3} | tail -n +2 | awk '{{print ${0} "_" ${1} "\t" ${2} "\t" ${5} }}' | sed s/chr//g | sort -k1,1 > {4}_tmp'''.format(chr_col, snp_col, pvalue_col, file, out, trait_col), shell=True)

    # Write header to final file
    subprocess.check_call('echo "chr\tpos\trsid\tpval\ttrait" > {0}'.format(out), shell=True)

    # Write rsids to final file
    subprocess.check_call('''join /users/mgloud/projects/gwas/data/rsid_mapping_sorted.txt {0}_tmp | grep rs | grep -v \; | grep -v , | sed s/\ /\\\t/g | sed s/_/\\\t/g | sort -k1,1n -k2,2n >> {0}'''.format(out), shell=True)

    subprocess.check_call("rm {0}_tmp".format(out), shell=True)


