echo -e "rsid\talt\tfreq\talt_freq\tFreqSE\tMinFreq\tMaxFreq\tbeta\tse\tpvalue\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tN\tID\tchr\tsnp_pos\tREF2\tALT2\tN_studies\tavsnp144\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tgwasCatalog" \
        > /users/mgloud/projects/gwas/data/prepared/MI_adjBMI_GENESISGuardian_Mixed.prepared.txt

tail -n +2 /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/MI_adjBMI_GENESISGuardian_Mixed_AllSNPs.txt | sort -k18,18n -k19,19n >> /users/mgloud/projects/gwas/data/prepared/MI_adjBMI_GENESISGuardian_Mixed.prepared.txt
bgzip -f /users/mgloud/projects/gwas/data/prepared/MI_adjBMI_GENESISGuardian_Mixed.prepared.txt
tabix -f -s 18 -b 19 -e 19 -S 1 /users/mgloud/projects/gwas/data/prepared/MI_adjBMI_GENESISGuardian_Mixed.prepared.txt.gz
