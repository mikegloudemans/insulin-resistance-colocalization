# Insulin Resistance Colocalization Analysis

Analysis performed by Mike Gloudemans

With contributions from Brunilda Balliu, Ivan Cárcamo-Oribe, Abhiram Rao, Joshua Knowles, Erik Ingelsson,
Matthew Durrant, Stephen Montgomery, and Thomas Quertermous

## Summary

This project contains the scripts required to perform the colocalization analyses described in the paper.

The top-level project directory should contain the folders `bin`, `data`, `output`, `tmp`, 
and `scripts`. All scripts should be
run from this top-level directory, or they'll be unable to locate the required files.

The complete analysis can be run sequentially from `pipe.sh` in the top-level project directory.
NOTE: The scripts can ONLY be run from this directory, as the file paths are specified relative to this
top level.


## Required Tools / Dependencies

The following analyses were performed using some other publicly available tools. To fully complete the
analyses, you will have to install or link these tools within the `bin` subfolder, or modify the
`pipe.sh` script to include the paths to the directories where these tools are installed.

* Tools for downloading and munging publicly available GWAS summary statistics
  are available at https://github.com/mikegloudemans/gwas-download/.
* Graphical visualization of colocalizations was performed using [LocusCompare](http://locuscompare.com) 
  (Liu et al. 2019)
* The analysis performed in this paper uses an integration of the publicly available tools 
  [FINEMAP](http://www.christianbenner.com/) (Benner et al. 2016)
  and [eCAVIAR](http://zarlab.cs.ucla.edu/tag/ecaviar/) (Hormozdiari et al. 2016), 
  available in a basic form at https://bitbucket.org/mgloud/production_coloc_pipeline/src. 

  An extended and hopefully easier-to-use pipeline with a greater variety of options and analyses
  will soon be available at https://github.com/mikegloudemans/ensemble_coloc_pipeline.

### Getting started

* An hg38-formatted version of the 
  [1000 Genomes VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/)
  is required for computing allele frequencies in a reference population.
* GWAS summary statistics are publicly available; consistently-formatted versions of these and other GWAS 
  can be [downloaded directly](https://github.com/mikegloudemans/gwas-download).
* GTEx v8 QTL association statistics can be downloaded from the [GTEx Portal](https://gtexportal.org/home/datasets). Some minor pre-processing will be required to run these scripts; this process is described [here](https://bitbucket.org/mgloud/production_coloc_pipeline/src).

### All required data

If you're having trouble accessing any of these data, please contact me (see Contact section
below) and I can quickly point you to the right location to obtain them, or set up a direct transfer.

To run all the scripts listed here, the `data` folder must contain all of the following files:

* `data/1KG`: 1KG VCF for hg38, publicly available for download as described above.
* `data/eqtls`: GTEx eQTLs for v8, and any other QTLs of interest, downloadable from GTEx Portal 
  as described above.
* `data/gwas`: All publicly available GWAS summary statistics, downloadable as described above. 
  If already re-formatted for colocalization, place them in a `formatted` subdirectory; if not, place 
  them in a `raw` subdirectory.

The rest of these files are not explicitly necessary, but were required to obtain the full results of the 
post-processing steps.

* `data/cadd`: Should contain a CADD file with VEP consequence predictions for every possible variant. May be 
  skipped if necessary.
* `data/hgnc`: The file `mart_export.txt` should contain a basic mapping of Ensembl gene IDs to HGNC names, 
  obtained through Biomart. This step can be omitted if such a file is not present.
* `data/ld`: Pairwise LD scores from 1K Genomes Phase 1, used for selecting LD buddies in the Variant Effect 
  Predictor annotation step, which can also be skipped if necessary.

## About the scripts

I've broadly organized the scripts into "pre-coloc", "colocalization", and "post-coloc" sections. Most of the
code in this project is focused on the "post-coloc" analysis, since the other two parts of the analysis
lean heavily on code from other projects (described and linked above).

Scripts to generate figures for the paper are in a dedicated folder.

Many of these scripts have wide-ranging applications but are currently geared towards our IR-specific
application. Feel free to contact me if you're trying to figure out how to gear a particular script
for your own analysis; some customization may be necessary but it's certainly doable.

I'm also working on a general Snakemake pipeline "post-coloc-toolkit", which I hope to make available
soon, and will work not only on our IR data but on any set of GWAS and QTL summary statistics! I will
link it from here when it's public.

## Contact

For any questions about this pipeline or about the colocalization-related analyses for this project, 
please contact Mike Gloudemans (mgloud@stanford.edu). I'll be glad to help you get these analyses up 
and running!

### Quick note on research ethics

We're making this code and the associated data available in hopes that it will help to further biomedical research.
Please be mindful of the ethical implications of your intended application, and use it for good :)
