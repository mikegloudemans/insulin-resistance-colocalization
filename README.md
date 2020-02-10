# Insulin Resistance Colocalization Analysis

##### Table of Contents  
[Headers](#headers)  
[Emphasis](#emphasis)

Analysis performed by Mike Gloudemans

## Summary

This project contains the scripts required to perform the colocalization analyses described in the paper.

The top-level project directory should contain the folders `bin`, `data`, `output`, `tmp`, and `scripts`. All scripts should be
run from this top-level directory, or they'll be unable to locate the required files.

The complete analysis can be run sequentially from `pipe.sh` in the top-level project directory.
NOTE: The scripts can ONLY be run from this directory, as the file paths are specified relative to this
top level.

## Methods (`scripts`)

Here I describe the individual steps in the analysis, all of which can be found in the `scripts` folder.

For the sake of clarity, I describe the individual steps in the order in
which they appear in `pipe.sh`, rather than alphabetizing them by filename.

### `scripts/pre_coloc`

Scripts to preprocess the GWAS and eQTL data, and to identify which loci will be tested for colocalization.

#### `scripts/pre_coloc/format_gwas/munge.sh`

Munge the original data files for IR-related traits into a format required to run the colocalization pipeline.
This step might not be necessary if the data are already downloaded in a munged format, as specified below in the
_Required Data_ section.

### `scripts/colocalization`

#### `scripts/colocalization/test_IR_colocs.py`

Just a wrapper script to run the main colocalization pipeline, which is available as a separate Github module as described
below.

### `scripts/post_coloc`

Various scripts needed to process the output of the colocalization algorithm so that
we can draw conclusions.

I eventually hope to incorporate these steps into a more user-friendly and generalized pipeline that can be applied
to an arbitrary set of GWAS.

#### `scripts/post_coloc/config`

Contains config files specifying the exact steps to be taken and filters to be used for the post-processing steps.

Also, specifies the location of the output directory where post-processed results should be placed.

#### `scripts/post_coloc/summarize_pre_test_statistics.R`

Counts the number of SNPs, loci, and genes available at each stage of filtering in the pipeline, and writes 
them to an output file located at `{output_directory}/full_coloc_results_qced.txt`.

Also, checks to make sure all tests were run as planned (i.e. the pipeline didn't exit early, or run into 
a formatting issue that forced it to skip over all instances of a certain GWAS). You'll need to examine
the output of the script to ensure that not a large number of tests were dropped, since this step
doesn't throw any actual error even if the results are unusual.

#### `scripts/post_coloc/aggregate_coloc_results.R`

Adds a few additional annotations, such as HGNC gene ID and rsid. Output will be placed at
`{output_directory}/clpp_results_{analysis_date}.txt`.

#### `scripts/post_coloc/last_qc_check.R`

The final QC check for any final errors in the file. You'll have to inspect the output yourself
to verify correctness. Specifically, we look for:
* Were the vast majority of candidate SNPs actually tested for colocalization? (Some dropout is normal due
  to minor eQTL / GWAS / reference genome discrepancies, as well as due to loci for which the eQTL gene lies
  too far from the GWAS locus to be confidently tested. But if more than ~5% of loci were dropped, this usually
  means something went wrong, and you'll have to look at the output from the colocalization pipeline itself to understand
  the error.)
* Does the final range of p-values for best GWAS SNPs per locus match our specified thresholds?
* Does the final range of p-values for best eQTL SNPs per locus match our specified thresholds?
* Does every individual GWAS and every individual eQTL study have a reasonable number of tested loci? (e.g. not 0)
* Are there a large number of loci for which not many SNPs (e.g. < 10) were present in the GWAS / eQTL overlap? (This would
  imply that a colocalization test is inappropriate for these loci.) If so, consult the accompanying table to see which files have many
  loci without many SNPs present.

#### `scripts/post_coloc/classify_results`

After the colocalization results have been aggregated into a single file and QC checks have been performed,
we then add an additional classification step for prioritizing loci with certain characteristics or from
specific high-priority GWAS, as described in the methods section of the paper. 

This analysis is performed using the R script in `scripts/post_coloc/classify_results/classify_results.R` and requires an
additional configuration file that specifies the exact classifications to be made. This config file (of which several examples
are already included) can be modified to 
work with virtually any GWAS and QTL files for arbitrary classification schema, but a description of that is outside of the
scope of this project. I plan to release that independently of this project and will link it here when published.

#### `scripts/post_coloc/get_variant_veps.R`

Fetch the Variant Effect Predictor annotations for all SNPs overlapping the loci of interest, as well as all "LD buddies" of
the lead GWAS SNPs (r2 < 0.8). This may not be necessary for all applications, but it gives us an initial sense of how many of the
candidate colocalizations are likely intergenic vs. how many overlap genic features such as coding regions, promoters, and/or introns.
Output file will be placed in `{output_directory}/vep_summary_{analysis_date}.txt`.

### `scripts/auxiliary`

A few other miscellaneous scripts for analyses that may be of interest but were not explicitly described in the paper.

#### `scripts/auxiliary/get_quantiles.R`

Empirically determine what levels of the modified CLPP score represent top 20% (weak) and 
top 5% (strong) colocalization in the list of candidate loci tested. This is done using the
CLPP results located at `{output_directory}/full_coloc_results_qced.txt`.


## Required Tools / Dependencies

The following analyses were performed using some other publicly available tools. To fully complete the
analyses, you will have to install or link these tools within the `bin` subfolder, or modify the
`pipe.sh` script to include the paths to the directories where these tools are installed.

* Tools for downloading and munging publicly available GWAS summary statistics
are available at https://github.com/mikegloudemans/gwas-download/.
* Graphical visualization of colocalizations was performed using [LocusCompare](http://locuscompare.com) (Liu et al. 2019)
* The analysis performed in this paper uses an integration of the publicly available tools [FINEMAP](http://www.christianbenner.com/) (Benner et al. 2016)
and [eCAVIAR](http://zarlab.cs.ucla.edu/tag/ecaviar/) (Hormozdiari et al. 2016), 
available in a basic form at https://bitbucket.org/mgloud/production_coloc_pipeline/src. 
An extended and hopefully easier-to-use pipeline with a greater variety of options and analyses, closer to what was used for this
paper, will soon be available at https://github.com/mikegloudemans/ensemble_coloc_pipeline.

## Required Data

This project makes uses of a variety of data tables, some already publicly available
and some internally generated. I'm currently exploring ways to just link all or most of the required
data files here, as well as key intermediate files and results, all as a single download. 
Until then, please contact me directly (see _Contact_ section below) and I'll share 
the relevant files directly, ASAP.

90% of the analysis, including the core colocalization analysis, can be completed using just the following
data files:

### Getting started

* An hg38-formatted version of the [1000 Genomes VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/)
is required for computing allele frequencies in a reference population.
* GWAS summary statistics are publicly available; consistently-formatted versions of these and other GWAS can be [downloaded directly](https://github.com/mikegloudemans/gwas-download).
* GTEx v8 QTL association statistics can be downloaded from the [GTEx Portal](https://gtexportal.org/home/datasets). Some minor pre-processing will be required to run these scripts; this process is described [here](https://bitbucket.org/mgloud/production_coloc_pipeline/src).

### All required data

To run all the scripts listed here, the `data` folder must contain all of the following files:

* `data/1KG`: 1KG VCF for hg38, publicly available for download as described above.
* `data/eqtls`: GTEx eQTLs for v8, and any other QTLs of interest, downloadable from GTEx Portal as described above.
* `data/gwas`: All publicly available GWAS summary statistics, downloadable as described above. If already re-formatted for colocalization, place them in a `formatted` subdirectory; if not, place them in a `raw` subdirectory.

The rest of these files are not explicitly necessary, but were required to obtain the full results of the post-processing steps.

* `data/cadd`: Should contain a CADD file with VEP consequence predictions for every possible variant. May be skipped if necessary.
* `data/hgnc`: The file `mart_export.txt` should contain a basic mapping of Ensembl gene IDs to HGNC names, obtained through Biomart. This step
can be skipped if such a file is not yet available.
* `data/ld`: Pairwise LD scores from 1K Genomes Phase 1, used for selecting LD buddies in the Variant Effect Predictor annotation step, which
can be skipped if necessary.

## Contact

For any questions about this pipeline or about the colocalization-related analyses for this project, 
please contact Mike Gloudemans (mgloud@stanford.edu). I'll be glad to help you get these analyses up and running!

## Quick note on research ethics

We're making this code and the associated data available in hopes that it will help to further biomedical research.
Please be mindful of the ethical implications of your intended application, and use it for good :)
