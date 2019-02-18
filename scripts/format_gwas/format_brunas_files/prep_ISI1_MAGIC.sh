echo -e "rsid\tchr\tsnp_pos\talt\tref\tmaf\tbeta\tse\tpvalue" \
	> /users/mgloud/projects/gwas/data/prepared/MAGIC_ISI_Model_1_AgeSexOnly.prepared.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/MAGIC_ISI_Model_1_AgeSexOnly.txt | cut -f2,3 --complement | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/gwas/data/prepared/MAGIC_ISI_Model_1_AgeSexOnly.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/MAGIC_ISI_Model_1_AgeSexOnly.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/MAGIC_ISI_Model_1_AgeSexOnly.prepared.txt.gz

