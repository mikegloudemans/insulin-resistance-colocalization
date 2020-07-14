## Generate colocalization heatmaps

The script `post_coloc/make_heatmaps.R` produces a heatmap from post-processed 
colocalization results.  either straight out of the script `post_coloc/aggregate_coloc_results.R`, 
or after further categorization by the script `post_coloc/classify_results/classify_results.R`.

In its simplest, default form, the script produces a heatmap in which the rows correspond to 
genes (and their genomic loci) and the columns correspond to combinations of QTL tissues with 
GWAS traits of interest. The colors in the heatmap indicate when colocalization occurred 
with eQTLs (blue), sQTLs (red), or both (purple), and the opacity of these colors indicates 
the relative evidence in support of colocalization at the [locus]-[gene]-[GWAS trait]-[QTL tissue] 
combination indicated by the cell’s position in the matrix.

Several easy customizations are available to facilitate visual inspection with heatmaps. 
One set of custom settings enables aggregating data across genes, GWAS traits, and tissues. 
Another allows subsetting the colocalization results into separate heatmaps based on a field of 
the results table (for example, a category created using the script `post_coloc/classify_results/classify_results.R`).  
The full set of options is detailed below in the "Parameters and settings" section.

### Parameters and settings

Parameters must be specified in a JSON-formatted config file, and its location passed to the 
script `post_coloc/make_heatmaps.R` as a command-line argument. The JSON file should contain 
an object (analogous to a Python dictionary) with parameter settings as described below. See 
the included example config file for an idea of how this file is created.

#### Required parameters

`"coloc_file"`: The path, relative to the current working directory, of the (TSV) file 
containing the colocalization results as generated using `post_coloc/aggregate_coloc_results.R`
or `post_coloc/classify_results/classify_results.R`.

`"out_folder_prefix"`: The path, relative to the current working directory, of the directory 
in which output heatmaps should be placed.

`"qtl_types"`: An object containing one or both of the keys `"eqtl"` and `"sqtl"`. Each 
of these keys should point to another nested object. The keys of that object must correspond 
to the values specified in the `"eqtl_file"` column of the `"coloc_file"` results file. 
The values paired with these keys should be the desired display names of this tissue 
within the output plot. NOTE: If an eQTL and sQTL file correspond to the same tissue, 
they must have the same display names specified, as shown in the example below.

Ex.
```
"qtl_types":
{
	"eqtl":
	{
		"data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz": "AdpSQ",
		"data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz": "AdpV",
		"data/eqtls/gtex_v8/Muscle_Skeletal.allpairs.txt.gz.eQTLs.txt.gz": "MuSk",
		"data/eqtls/gtex_v8/Liver.allpairs.txt.gz.eQTLs.txt.gz": "Liv",
		"data/eqtls/gtex_v8/Pancreas.allpairs.txt.gz.eQTLs.txt.gz": "Panc"
	},
	"sqtl":
	{
		"data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz": "AdpSQ",
		"data/sqtls/gtex_v8/Liver.sQTLs.txt.gz": "Liv",
		"data/sqtls/gtex_v8/Muscle_Skeletal.sQTLs.txt.gz": "MuSk",
		"data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz": "AdpV",
		"data/sqtls/gtex_v8/Pancreas.sQTLs.txt.gz": "Panc"
	}
},
```

`"tissue_order"`:  An array (analogous to a Python list) containing the order in which the 
tissue names should be displayed in the output plot(s).

Ex.

```
"tissue_order": ["AdpSQ", "AdpV",  "MuSk", "Liv", "Panc"],
```

`"gwas_labels"`: An object, mapping the GWAS traits in the `base_gwas_file` column of the `"coloc_file"` 
results to the desired display names of these traits within the output plot. Note that 
multiple GWAS files can be collapsed into the same label; in the example below, four different 
T2D files have been collapsed into a single T2D GWAS category for display in the output heatmap.

Ex.
```
"gwas_labels":
{
	"data/gwas/formatted/sumstats/hg38/HDL_GLGC_Expanded/HDL_GLGC_Expanded.txt.gz": "HDL",
	"data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Xue_2018/Type-2-Diabetes_Xue_2018.txt.gz": "T2D",
	"data/gwas/formatted/sumstats/hg38/BMI_GIANT_2018/BMI_GIANT_2018.txt.gz": "BMI",
	"data/gwas/formatted/sumstats/hg38/TG_GLGC_Expanded/TG_GLGC_Expanded.txt.gz": "TG",
	"data/gwas/formatted/sumstats/hg38/WHR-adj-BMI_GIANT_2018/WHR-adj-BMI_GIANT_2018.txt.gz": "WHR",
	"data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Spracklen_2020/Type-2-Diabetes_Spracklen_2020.txt.gz": "T2D",
	"data/gwas/formatted/sumstats/hg38/T2D_Mahajan_Europeans/T2D_Mahajan_Europeans.txt.gz": "T2D",
	"data/gwas/formatted/sumstats/hg38/Type-2-Diabetes_Suzuki_2019/Type-2-Diabetes_Suzuki_2019.txt.gz": "T2D ",
	"data/gwas/formatted/sumstats/hg38/MI_adjBMI_European/MI_adjBMI_European.txt.gz": "ISI_Genesis",
	"data/gwas/formatted/sumstats/hg38/MAGIC_ISI_Model_2_AgeSexBMI.txt/MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz": "ISI_MAGIC",
	"data/gwas/formatted/sumstats/hg38/FastGlu_MAGIC_Europeans/FastGlu_MAGIC_Europeans.txt.gz": "FastGlu",
	"data/gwas/formatted/sumstats/hg38/FastInsu_adjBMI_MAGIC_Europeans/FastInsu_adjBMI_MAGIC_Europeans.txt.gz": "FastInsu"
},
```

`"gwas_order"`: An array containing the order in which the GWAS trait names should be 
displayed in the output plot(s).
```
"gwas_order": ["ISI_Genesis", "ISI_MAGIC", "FastInsu", "FastGlu", "T2D", "WHR","HDL","TG","BMI"],
```

`"file_strata"`: A list of JSON objects. For each object in this list, an additional 
set of heatmaps will be generated. Each JSON object should contain the following fields, 
which instruct the script on which files to create and how to create them:

* `"out_dir"` (required): a suffix to add to the output folder prefix specified 
by `"out_folder_prefix"` at the top level JSON object.
* `"split_factors"` (required): A list of columns on which to stratify the output 
files. If more than one column is specified, an output heatmap will be created 
for every possible combination of factors in these columns.
* `"concise"` (required): If `"True"`, remove all rows from the output heatmap 
in which no colocalization was found. Otherwise, include all rows in which at
 least one test was performed, even if no colocalization was detected.

Ex.

```    
"file_strata" : [         
{
                        "out_dir": "concise/non-pancreas-stratified/t2d_collapsed/step3/",
                        "split_factors": ["step3"],
                        "concise": "True",
                        "gwas_remapping":
                        {
                                "remaps":
                                {
                                        "T2D-Mahajan": "T2D",
                                        "T2D-Suzuki": "T2D",
                                        "T2D-Xue": "T2D",
                                        "T2D-Spracklen": "T2D"
                                },
                                "new_order": ["ISI_Genesis", "ISI_MAGIC", "FastInsu", "FastGlu", "T2D","WHR","HDL","TG","BMI"]
                        }
                },
                {
                        "out_dir": "concise/pancreas-stratified/step1/",
                        "split_factors": ["step1", "pancreas"],
                        "concise": "True"
                }
]
```

#### Optional parameters

`"chunk_size"`: The maximum number of rows to be displayed in a single heatmap. 
If the number of rows exceeds this number, the heatmap will be broken into chunks 
of this size and output as separate files. If not specified, the default value is 100.

`"mark_dubious_results"`: If `"True"`, display an "X" or another mark in cells in the 
heatmap containing results for which the GWAS or eQTL pvalue did not pass some specified 
cutoff, regardless of the colocalization strength.

`"cluster"`: If `"True"`, cluster rows of the output heatmap using the Jaccard 
index similarity metric. If "False" or unspecified, just order them by genomic 
location. Columns remain in the same order regardless.

`"y_axis_collapse"`: 

* If `"genes"`, collapse rows for all genes tested at each locus into a single row, 
indicating whether any colocalization occurred at that locus in the given GWAS-tissue combo. 

* If anything else or if not specified, maintain a row for every gene tested at the 
locus in any GWAS-tissue combo.

`"x_axis_collapse"`:
 
* If `"gwas"`, collapse columns for all GWAS traits tested in a tissue, indicating 
whether any colocalization occurred in this tissue at a given locus-gene combo. 
The columns will then be just one for each QTL tissue.

* If `"tissues"`, collapse columns for all tissues tested with each GWAS trait, 
indicating whether any colocalization occurred in this GWAS at a given locus-gene 
combo. The columns will then be just one for each GWAS.

* If `"tissues-gwas"`, collapse columns for all GWAS traits tested in all tissues, 
indicating whether any colocalization occurred at a given locus-gene combo. 
There will be just one column if this option is chosen.

* If anything else or if not specified, there will be a column for every GWAS-tissue combination.

If you want to run this code many times to create heatmaps with different settings, we 
recommend making a template JSON config script, and then manipulating it programmatically to 
run the script multiple times with the desired settings. 

