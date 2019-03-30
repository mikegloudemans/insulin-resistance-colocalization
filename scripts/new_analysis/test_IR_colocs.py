import operator
import subprocess
import json
import sys
import time

def main():
    # Reset things fresh on each run, so we're not mixing results
    subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/insulin-full-analysis-rerun/*", shell=True)

    kept_data = set([])
    with open("/users/mgloud/projects/insulin_resistance/output/snps_to_test.txt") as f:
        all_data = []
        f.readline()
        for line in f:
            data = line.strip().split()
            kept_data.add((data[0], data[1], data[7]))
    kept_data = list(kept_data)

    kept_data = sorted(kept_data, key=operator.itemgetter(2))

    # Then for every locus in the "kept data"...
    for i in range(len(kept_data)):

        test = kept_data[i]
        print test
        
        temp = json.loads(template)
        temp["snp_list_file"] = "/users/mgloud/projects/insulin_resistance/tmp/snp_list{0}.txt".format(i)

        # Add locus to SNP list...but only once for each gene
        with open("/users/mgloud/projects/insulin_resistance/tmp/snp_list{0}.txt".format(i), "w") as w:
            w.write("{0}\t{1}\t{2}\n".format(test[0], test[1], test[2]))
               
        # Write config file to the appropriate directory
        with open("/users/mgloud/projects/insulin_resistance/tmp/ir_config{0}.config".format(i), "w") as w:
            json.dump(temp, w)

        # Run the test
        subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/insulin_resistance/tmp/ir_config{0}.config 1 &".format(i), shell=True)

        while int(subprocess.check_output('''ps -ef | grep "python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/insulin_resistance/tmp/ir_config" | wc -l''', shell=True)) > 8:
            time.sleep(5)

template = '''
{
        "out_dir_group": "insulin-full-analysis-rerun",

       "gwas_experiments": 
	{
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_LDL_Mixed.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_MI_adjBMI.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_MI.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_T2D_DIAGRAM_European.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/bonus_sumstats/GWAS_WHR_European.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_2hrGlu_MAGIC_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_BMI_GIANT_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_CHD_CARDIoGRAMplusC4D_Mixed.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_FastGlu_MAGIC_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_FastInsu_adjBMI_MAGIC_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_FastInsu_MAGIC_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_HDL_GCLC_Mixed.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_TG_GCLC_Mixed.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_WHRadjBMI_GIANT_Europeans.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"},
            "/users/mgloud/projects/insulin_resistance/data/gwas/formatted/brunas_sumstats/GWAS_MAGIC_ISI_Model_1_AgeSexOnly.txt.txt.gz": {"ref": "1kgenomes", "gwas_format": "pval_only"}

	},
	
	"eqtl_experiments":	
	{
            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/Adipose_Subcutaneous.allpairs.txt.gz": {"ref": "1kgenomes", "eqtl_format": "effect_size"},
            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/Adipose_Visceral_Omentum.allpairs.txt.gz": {"ref": "1kgenomes", "eqtl_format": "effect_size"},
            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/Liver.allpairs.txt.gz": {"ref": "1kgenomes", "eqtl_format": "effect_size"},
            "/users/mgloud/projects/brain_gwas/data/eqtls/gtex_v7/Muscle_Skeletal.allpairs.txt.gz": {"ref": "1kgenomes", "eqtl_format": "effect_size"}
	},

	"eqtl_threshold": 
		1,

	"selection_basis": 
		"snps_from_list",

	"snp_list_file":
                "/users/mgloud/projects/insulin_resistance/tmp/snp_list.txt",

	"methods": 
	{
		"finemap":{}
	},

        "ref_genomes": 
	{
		"1kgenomes": 
		{
			"file": 
				"/mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",

                	"af_attribute": 
				"AF",

                        "N": 
				2504
	        }
        }
}
'''

if __name__ == "__main__":
    main()
