require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

# Format colocalization results table

# For this script, we assume that the accompanying script `summarize_pre_test_statistics.R`
# has already been run for the colocalization run of interest.

# These results have been filtered to include only the relevant GWAS traits and to remove
# any colocalization runs that didn't pass significance after intersection with the eQTL
# summary statistics.

### Load filtered and QC'ed COLOC results from the previous script.

config_file = commandArgs(trailing=TRUE)[1]
config = fromJSON(file=config_file)

results = read.table(paste(config$out_dir, "full_coloc_results_qced.txt", sep="/"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
results$chr = sapply(as.character(results$ref_snp), function(x) {strsplit(x, "_")[[1]][1]})
results$pos = sapply(as.character(results$ref_snp), function(x) {strsplit(x, "_")[[1]][2]})

# We collect and add in the necessary metadata to fill out the specified
# data output format.

# Get HGNC gene names
genes = read.table("data/hgnc/mart_export.txt", header=TRUE, sep="\t")
genes = genes[,-2]
names(genes) = c("ensembl", "hgnc")
genes = genes[!duplicated(genes$ensembl),]
results = merge(genes, results, all.y=TRUE)

results$hgnc = as.character(results$hgnc)
results$ensembl = as.character(results$ensembl)

results$hgnc[is.na(results$hgnc)]=results$ensembl[is.na(results$hgnc)]
results$hgnc[results$hgnc==""]=results$ensembl[results$hgnc==""]


# Load rsIDs -- here, we'll use a simple lookup in the 1Kgenomes file
rsids = results[!duplicated(results$ref_snp),][c("ref_snp", "chr", "pos")]
rsids$rsid = sapply(1:dim(rsids)[1], function(i)
	{
		chr = rsids$chr[i]
		pos = rsids$pos[i]
		query = system(paste0("tabix data/1KG/ALL.chr", chr, "_GRCh38.genotypes.20170504.vcf.gz ", chr, ":", pos, "-", pos), intern=TRUE)
		if (length(query) > 1)
		{
			vars = sapply(1:length(query), function(i) {strsplit(query[[i]], "\\t")[[1]][3]})
			vars2 = vars[grepl("rs", vars)]
			if (length(vars2) > 0)
			{
				return(paste(vars2, collapse="_"))
			}
			return(paste(vars, collapse="_"))
		}
		if (length(query) == 0)
		{
			print("No variant found at this locus. Please check the code for errors.")
			return(rsids$ref_snp[i])
		}
		rsid = strsplit(query, "\\t")[[1]][3]
		return(rsid)
	}
)
results = merge(rsids[c("ref_snp", "rsid")], results, by="ref_snp")

results$min_gwas_pval = 10^-(results$neg_log_gwas_pval)
results$min_eqtl_pval = 10^-(results$neg_log_eqtl_pval)

### Output our results file for further in-depth analysis

results_output = results[c("rsid", "chr", "pos", "locus", "min_gwas_pval", "min_eqtl_pval", "gwas_short", "eqtl_short", "gwas_trait", "eqtl_file", "feature", "ensembl", "hgnc", "clpp", "clpp_mod", "n_snps", "all_sig_gwas", "all_sig_eqtl", "base_gwas_file")]
write.table(results_output, paste0(config$out_dir, "/clpp_results_", config$analysis_date, ".txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
