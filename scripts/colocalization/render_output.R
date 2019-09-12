# Wrapper script to render the Rmarkdown scripts for a given config file

args = commandArgs(trailingOnly=TRUE)
config_file = args[1]

require(rmarkdown)

render("/users/mgloud/projects/insulin_resistance/scripts/colocalization/summarize_pre_test_statistics.Rmd")
render("/users/mgloud/projects/insulin_resistance/scripts/colocalization/aggregate_coloc_results.Rmd")
