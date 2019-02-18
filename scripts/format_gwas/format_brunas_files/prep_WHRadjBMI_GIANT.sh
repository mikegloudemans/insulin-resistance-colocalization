# WARNING! We're not yet sure if the ref/alt direction is correct for these GWAS files.

echo -e "rsid\tchr\tsnp_pos\tref\talt\tref_freq\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/WHRadjBMI_GIANT_Europeans_AllSNPs.txt | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Europeans_AllSNPs.prepared.txt.gz


echo -e "rsid\tchr\tsnp_pos\tref\talt\tref_freq\tbeta\tse\tpvalue\tN" \
> /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/WHRadjBMI_GIANT_Mixed_AllSNPs.txt | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/WHRadjBMI_GIANT_Mixed_AllSNPs.prepared.txt.gz


