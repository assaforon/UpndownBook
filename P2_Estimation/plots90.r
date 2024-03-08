cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(ggplot2)
theme_set(theme_bw(16)) # To be changed to bw later
library(patchwork)

load(file.path(outdir, 'grandsim90w.RData'))	
e90w = ls(pat='est[bk]w9')
#e90w = e90w[!grepl('[57]$', e90w) ]

#### Constants

wid = 12
hgt = 7

desnames = c('BCD', "K-row")
colors4 = c('grey65', 'grey40', 'grey80', 'black')
inames = c('Coverage', 'Width')

gpoints = c("dm48", "dyna", "cir",  "ir",   "dual")

# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')

p90stack = combo(e90w, finites = FALSE)
p90stack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]
i90stack = combo(e90w, atomfun = imetrix, outnames = inames)
i90stack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]

p90stack[ , baseframe := gsub('[_]slow', '', Framework) ]
p90mash = merge(p90stack[Metric == 'MAE90' & estimate=='cir' & !grepl('slow', Framework)],
				p90stack[Metric == 'MAE90' & estimate=='cir' & grepl('slow', Framework)],
					all = FALSE, by = c('Design', 'Metric', 'baseframe') )
					
i90stack[ , baseframe := gsub('[_]slow', '', Framework) ]
i90mash = merge(i90stack[Metric == 'Coverage' & estimate=='cir' & !grepl('slow', Framework)],
		i90stack[Metric == 'Coverage' & estimate=='cir' & grepl('slow', Framework)],
					all = FALSE, by = c('Design', 'Metric', 'baseframe') ) 

	y1max = max(c(p90mash$Value.x, p90mash$Value.y))

labelz = scale_x_continuous(labels = c(paste('BCD,',c('High', 'Low'),'Start'), paste('K-row,',c('High', 'Low'),'Start') ) )

rspace = theme(plot.margin = unit(c(0,50,0,0), "pt"))	
	p90mash[ , align := 1:4]	
	p1 <- ggplot(p90mash, aes(x=align, y=Value.y)) + geom_point(color='grey', size = 5) +
			geom_segment(aes(xend=align, yend=Value.x), arrow = arrow(length=unit(0.03,'npc')) , lineend='round',linejoin='round', linewidth=1.5) + ylim(0,y1max) + geom_point(aes(y=Value.x), size = 5) + labs(x='', y='MAE90') + labelz
	
	i90mash[ , align := 1:4]	
	p2 <- ggplot(i90mash, aes(x=align, y=100*Value.y)) + geom_point(color='grey', size = 5) +
			geom_segment(aes(xend=align, yend=100*Value.x), arrow = arrow(length=unit(0.03,'npc')) , lineend='round',linejoin='round', linewidth=1.5) + geom_point(aes(y=100*Value.x), size = 5) + labs(x='', y='Interval Coverage (%)') + labelz +
			geom_hline(yintercept = 90) + geom_hline(yintercept = 85, lty=2)
	
ggsave(p1+rspace+p2+rspace, file = file.path(outdir, 'sim_quickstart.pdf'),
				width = 13, height=6)



	stop('a')
	
	ggplot(p90mash, aes(Value.y, Value.x, color = Design)) + geom_point(size=5) +
			geom_abline(intercept=0, slope=1) + xlim(0.03,y1max) +  ylim(0.03,y1max) +
			scale_color_manual (values = colors4[c(1,4)]) + 
			labs(x='MAE90, Naive Design', y='MAE90, with Quick-Start')

	y2range = 100 * range(c(i90mash$Value.x, i90mash$Value.y))				

	p2 <- ggplot(i90mash, aes(100*Value.y, 100*Value.x, color = Design)) + 
		geom_point(size=3) + geom_vline(xintercept=90) +
			geom_hline(yintercept=90) + xlim(y2range) + ylim(y2range) +
			scale_color_manual (values = colors4[c(1,4)]) + 
			labs(x='CI Coverage, Naive Design', y='CI Coverage, with Quick-Start', color = NULL)
			
			
					
					



