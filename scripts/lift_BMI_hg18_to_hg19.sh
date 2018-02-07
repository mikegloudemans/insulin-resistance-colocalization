# WARNING! We're not yet sure if the ref/alt direction is correct for these GWAS files.

echo -e "rsid\tchr\tsnp_pos\tgene\talt\tref\talt_freq\tbeta\tse\tvariance_explained\tN\tpvalue" \
	> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(tail -n +2 /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits.txt | sort -k1,1) | sort -k2,2 -k3,3n | cut -f4,5 --complement >> /users/mgloud/projects/insulin_resistance/data/BMI_GIANT_Europeans_Hits_hg19.txt

