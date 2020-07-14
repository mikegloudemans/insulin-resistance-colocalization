require(rjson)

out_dir = "output/post_coloc/plots/heatmaps/" 
#config_template = "scripts/post_coloc/heatmaps/heatmap_config_template.config"
config_template = "scripts/post_coloc/heatmaps/heatmap_config_template_test.config"
json_out_file = "scripts/post_coloc/heatmaps/tmp_config.config"

rscript_file = "scripts/post_coloc/heatmaps/make_heatmaps.R"

config = fromJSON(file=config_template)

for (cluster_mode in c("True", "False"))
{
	for (y_axis_collapse in c("genes", "none"))
	{
		for (x_axis_collapse in c("tissues", "gwas", "tissues_gwas", "none"))
		{
			config_dup = config
			config_dup$cluster = cluster_mode
			config_dup$y_axis_collapse = y_axis_collapse
			config_dup$x_axis_collapse = x_axis_collapse
			heatmap_out_base = sprintf("%s//yaxis_collapsed_%s/x_axis_collapsed_%s/cluster_%s/", out_dir, y_axis_collapse, x_axis_collapse, cluster_mode)	
			config_dup$out_folder_prefix = heatmap_out_base

			writeLines(toJSON(config_dup), json_out_file)
			print(sprintf("%s %s %s", cluster_mode, y_axis_collapse, x_axis_collapse))
			system(sprintf("Rscript %s %s", rscript_file, json_out_file))
		}
	}
}
