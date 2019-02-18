echo -e "chr\tsnp_pos\talt\tref\tbeta\tse\tpvalue\tN" \
	> /users/mgloud/projects/gwas/data/prepared/T2D_DIAGRAM_European.prepared.txt
tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/T2D_DIAGRAM_European_AllSNPs.txt | sed s/:/\\t/g | sort -k1,1n -k2,2n >> /users/mgloud/projects/gwas/data/prepared/T2D_DIAGRAM_European.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/T2D_DIAGRAM_European.prepared.txt
tabix -f -s 1 -b 2 -e 2 -S 1 /users/mgloud/projects/gwas/data/prepared/T2D_DIAGRAM_European.prepared.txt.gz
