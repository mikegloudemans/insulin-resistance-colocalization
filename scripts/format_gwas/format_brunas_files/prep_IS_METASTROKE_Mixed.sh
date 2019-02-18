echo -e "rsid\tchr\tsnp_pos\talt\tref\tFreq1\tFreqSE\tbeta\tse\tpvalue\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\told_chromosome\told_position" \
	> /users/mgloud/projects/gwas/data/prepared/IS_METASTROKE_Mixed.prepared.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/IS_METASTROKE_Mixed_AllSNPs.txt | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/gwas/data/prepared/IS_METASTROKE_Mixed.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/IS_METASTROKE_Mixed.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/IS_METASTROKE_Mixed.prepared.txt.gz

