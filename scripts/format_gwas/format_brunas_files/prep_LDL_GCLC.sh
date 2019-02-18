echo -e "rsid\tchr\tsnp_pos\talt\tref\tbeta\tse\tN\tpvalue\talt_freq" \
	> /users/mgloud/projects/gwas/data/prepared/LDL_GCLC_Mixed.prepared.txt
join -1 3 -2 1 -t $'\t' <(zcat /users/mgloud/projects/gwas/data/sorted_hg19_snp150.txt.gz) <(zcat /mnt/lab_data/montgomery/shared/gwas/Blood-Lipids_Willer_2013/jointGwasMc_LDL.txt.gz | tail -n +2 | cut -f1,2 --complement | sort -k1,1) | sort -k2,2 -k3,3n >> /users/mgloud/projects/gwas/data/prepared/LDL_GCLC_Mixed.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/LDL_GCLC_Mixed.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/LDL_GCLC_Mixed.prepared.txt.gz

