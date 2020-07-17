suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(require(rjson)))
suppressWarnings(suppressMessages(require(gganimate)))

# Load color scheme
source("scripts/post_coloc/figures/color_scheme.R")

# First get two plot frames using the existing code

# I saved the shot, so just go from there for now

# save.image("animation_frames.RData")


frame1 = frame1[frame1$locus < 50,]
frame2 = frame2[frame2$locus < 50,]
frame1$frame = "frame1"
frame2$frame = "frame2"

all_frames = rbind(frame1, frame2)

all_frames$animation_grouping = paste(all_frames$locus, all_frames$x_factor)

plot=ggplot(data = all_frames, mapping = aes(x = x_factor, y = y_factor)) +
	geom_tile(mapping = aes(fill=coloc_class, group=animation_grouping),color="black") +
        scale_fill_manual(name="CLPP", values = color_scheme, drop=FALSE) +
	theme(axis.title = element_blank(),
		axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5),
		axis.text.y = element_text(size=12),
		legend.position = 'top',
		legend.box = "vertical",
		legend.text = element_text(size=15),
		legend.title = element_text(size = 12)) +
		guides(color=guide_legend(title = "", nrow = 2)) +
        	scale_x_discrete(drop=FALSE)

anim = plot + transition_states(frame, transition_length = 3, state_length = 5) #+
	#view_follow()
animate = animate(anim, height = 800, width = 800)
anim_save("test_anim.gif")

################# Get to all this later

row_count = length(unique(frame1$y_factor))
max_chunk = floor(row_count / chunk_size + 1)

genes_per_locus = frame1 %>% group_by(locus) %>% summarize(genes_at_locus=length(unique(y_factor))) %>% arrange(locus)

num_cols = length(levels(frame1$x_factor))

if (("x_axis_collapse" %in% names(config)) && ((config$x_axis_collapse == "tissues") || (config$x_axis_collapse == "tissues-gwas")))
{
	num_tissues = 1
} else
{
	num_tissues = length(unique(coloc_res$tissue))
}

num_vert_bars = num_cols / num_tissues - 1
num_rows = length(unique(frame1$y_factor))

if ((("cluster" %in% names(config)) && (config$cluster == "True")) ||
	(("y_axis_collapse" %in% names(config)) && (config$y_axis_collapse == "genes")))
{
	# It doesn't make sense to separate loci if they're clustered
	# It also doesn't make sense to separate loci if each row is a distinct locus
	num_horz_bars = 0
} else
{
	num_horz_bars = length(unique(tmp_chunk$locus))
}

horz_breaks = cumsum(genes_per_locus$genes_at_locus)

margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$y_factor))), 0))

plot=plot_coloc_results_function(data = tmp_chunk)

if (num_vert_bars != 0)
{
	my.vertical.lines<-data.frame(x=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, y = rep(0.5, num_vert_bars),
				      xend=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, yend = rep(num_rows + 0.5, num_vert_bars))
	plot = plot + geom_segment(data=my.vertical.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
}
if (num_horz_bars != 0)
{
	my.horizontal.lines<-data.frame(x=rep(0.5, num_horz_bars), y=num_rows-horz_breaks+0.5,
					  xend=rep(num_cols+0.5, num_horz_bars), yend=num_rows - horz_breaks+0.5)
	plot = plot + geom_segment(data=my.horizontal.lines, aes(x,y,xend=xend, yend=yend), size=0.25, inherit.aes=F)
}
plot

if(save_plots)
{
	ggsave(filename = paste0(plot_out_dir, '/', out_sub_folder, '/CLPP_group_',levels(coloc_res[["split_column"]])[i],'.part', chunk, '.pdf'), plot = plot, width = margin_approx_size+(col_width*num_cols), height = margin_approx_size+(row_height*num_rows), limitsize = F)
}


						
