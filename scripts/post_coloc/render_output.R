# Wrapper script to render the Rmarkdown scripts for a given config file

args = commandArgs(trailingOnly=TRUE)
config_file = args[1]

require(rmarkdown)

render("scripts/post_coloc/summarize_pre_test_statistics.Rmd")
render("scripts/post_coloc/aggregate_coloc_results.Rmd")
