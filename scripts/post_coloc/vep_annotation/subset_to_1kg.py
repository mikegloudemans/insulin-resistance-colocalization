import subprocess

melted_vep_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_2020-05-11.txt"
out_vep_file = "output/post_coloc/2020-05-11/refiltered/eqtls_and_sqtls/vep_summary_melted_1kg_filtered_2020-05-11.txt"

with open(out_vep_file, "w") as w:
    with open(melted_vep_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            chrom = data[0]
            if "_" in chrom:
                continue
            pos = int(data[1])
            ref = data[4]
            alt = data[5]
            command = "tabix data/1KG/ALL.chr{0}_GRCh38.genotypes.20170504.vcf.gz {0}:{1}-{2}".format(chrom, pos, pos+1)

            output = subprocess.check_output(command, shell=True).strip().split("\n")
            for out in output:
                if out == "":
                    continue
                vcf_data = out.strip().split("\t")

                vcf_pos = int(vcf_data[1])
                vcf_ref = vcf_data[3]
                vcf_alt = vcf_data[4]

                if vcf_pos != pos:
                    continue

                # Only keep the variant if it's actually found in 1k genomes
                if (vcf_ref == ref and vcf_alt == alt) or \
                        (vcf_ref == alt and vcf_alt == ref):
                    w.write(line)


