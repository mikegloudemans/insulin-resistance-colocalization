import operator
import subprocess
import json
import sys
import time

def main():
    # Reset things fresh on each run, so we're not mixing results
    subprocess.call("rm -rf output/colocalization/raw/ir-v8/*", shell=True)
    subprocess.call("rm -rf tmp", shell=True)
    subprocess.call("mkdir tmp", shell=True)

    kept_data = []
    with open("output/test_snps/ir-v8/ir-v8_all-gwas_gtex-ir_gwas-pval1e-05_eqtl-pval1e-05_gwas-window500000_eqtl-window0_coloc-tests.txt") as f:
        all_data = []
        f.readline()
        for line in f:
            data = line.strip().split()
            kept_data.append(data)
    kept_data = list(kept_data)

    # Then for every locus in the "kept data"...
    for i in range(len(kept_data)):

        test = kept_data[i]
        print test
        
        temp = json.loads(template)
        temp["snp_list_file"] = "tmp/snp_list{0}.txt".format(i)

        # Add locus to SNP list...but only once for each gene
        with open("tmp/snp_list{0}.txt".format(i), "w") as w:
            w.write("{0}\t{1}\t{2}\n".format(test[0], test[1], test[7]))
               
        # Add corresponding gwas experiment to the list, if not already present
        temp["gwas_experiments"][test[2]] = {"ref": "1kgenomes", "gwas_format": "pval_only"}
        if test[2] != test[4]:
            temp["gwas_experiments"][test[2]]["traits"] = [test[4]]

        # Add corresponding eQTL tissue to the list
        temp["eqtl_experiments"][test[3]] = {"ref": "1kgenomes", "eqtl_format": "effect_size"}

        # Write config file to the appropriate directory
        with open("tmp/ir_config{0}.config".format(i), "w") as w:
            json.dump(temp, w)

        # Run the test
        subprocess.call("python bin/colocalization_pipeline/dispatch.py tmp/ir_config{0}.config 1 &".format(i), shell=True)

        while int(subprocess.check_output('''ps -ef | grep "python bin/colocalization_pipeline/dispatch.py tmp/ir_config" | wc -l''', shell=True)) > 20:
            time.sleep(1)

    # Clean up
    subprocess.call("rm -rf tmp", shell=True)
    subprocess.call("mkdir tmp", shell=True)

template = '''
{
        "out_dir": "output/colocalization/raw/ir-v8",

       "gwas_experiments": 
	{
	},
	
	"eqtl_experiments":	
	{
	},

	"eqtl_threshold": 
		1,

	"selection_basis": 
		"snps_from_list",

	"snp_list_file":
                "tmp/snp_list{0}.txt",

	"methods": 
	{
		"finemap":{}
	},

        "ref_genomes": 
	{
		"1kgenomes": 
		{
			"file": 
                                "data/1KG/ALL.chr{0}_GRCh38.genotypes.20170504.vcf.gz",

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
