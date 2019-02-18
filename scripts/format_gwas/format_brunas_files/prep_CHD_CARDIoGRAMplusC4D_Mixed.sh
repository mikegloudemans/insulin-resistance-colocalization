echo -e "rsid\tchr\tsnp_pos\talt\tref\talt_freq\tmedian_info\tmodel\tbeta\tse\tpvalue\thet_p\tN" \
	> /users/mgloud/projects/gwas/data/prepared/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt
tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/CHD_CARDIoGRAMplusC4D_Mixed_AllSNPs.txt | sort -k2,2n -k3,3n >> /users/mgloud/projects/gwas/data/prepared/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt
tabix -f -s 2 -b 3 -e 3 -S 1 /users/mgloud/projects/gwas/data/prepared/CHD_CARDIoGRAMplusC4D_Mixed.prepared.txt.gz
