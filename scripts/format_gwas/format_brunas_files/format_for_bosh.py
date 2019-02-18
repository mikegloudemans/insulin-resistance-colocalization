import glob
import subprocess

for file in glob.glob("/users/mgloud/projects/gwas/data/prepared/*.gz"):

    if "QTL" not in file:
        continue

    print "Reformatting", file
    base =  file.replace(".gz", "").replace(".prepared.", ".").split("/prepared/")
    out = base[0] + "/boshs_format/GWAS_" + base[1]

    subprocess.check_call("python remap_rsids.py {0} {1}".format(file, out), shell=True)

#subprocess.check_call(["bash", "name_boshs_files.sh"])
