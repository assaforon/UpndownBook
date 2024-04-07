cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(ggplot2)
theme_set(theme_bw(18)) 
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


#### Creating a regular vs. quick-start dataset - commented out b/c already saved.

# p90mash = merge(p90stack[Metric == 'MAE90' & estimate=='cir' & !grepl('slow', Framework)],
				# p90stack[Metric == 'MAE90' & estimate=='cir' & grepl('slow', Framework)],
					# all = FALSE, by = c('Design', 'Metric', 'baseframe') )
					
# i90stack[ , baseframe := gsub('[_]slow', '', Framework) ]
# i90mash = merge(i90stack[Metric == 'Coverage' & estimate=='cir' & !grepl('slow', Framework)],
		# i90stack[Metric == 'Coverage' & estimate=='cir' & grepl('slow', Framework)],
					# all = FALSE, by = c('Design', 'Metric', 'baseframe') ) 

# #------------------ Plotting quick-start
# ### plot constants

# y1max = max(c(p90mash$Value.x, p90mash$Value.y))

# labelz = scale_x_continuous(labels = c(paste('BCD,',c('High', 'Low'),'Start'), paste('K-row,',c('High', 'Low'),'Start') ) )

# rotate = theme(axis.text.x = element_text(angle = 45, hjust=1) )

# rspace = theme(plot.margin = unit(c(0,30,0,0), "pt"))	

# ### The actual plots

	# p90mash[ , align := 1:4]	
	# p1 <- ggplot(p90mash, aes(x=align, y=Value.y)) + geom_point(color='grey', size = 5) +
			# geom_segment(aes(xend=align, yend=Value.x), arrow = arrow(length=unit(0.03,'npc')) , lineend='round',linejoin='round', linewidth=1.5) + ylim(0,y1max) + geom_point(aes(y=Value.x), size = 5) + labs(x='', y='MAE90 (spacing units)') + labelz + rotate
	
	# i90mash[ , align := 1:4]	
	# p2 <- ggplot(i90mash, aes(x=align, y=100*Value.y)) + geom_point(color='grey', size = 5) +
			# geom_segment(aes(xend=align, yend=100*Value.x), arrow = arrow(length=unit(0.03,'npc')) , lineend='round',linejoin='round', linewidth=1.5) + geom_point(aes(y=100*Value.x), size = 5) + labs(x='', y='Interval Coverage (%)') + labelz +
			# geom_hline(yintercept = 90) + geom_hline(yintercept = 85, lty=2) + rotate
	
### Arranging and saving; the syntax is from the "patchwork" package

### Commenting out b/c it's already exported and uploaded; uncomment to re-run

#ggsave(p1+rspace+p2+rspace, file = file.path(outdir, 'sim_quickstart.pdf'),
#				width = 13, height=7)


#---------------------- "Standard" Plots

p90stack = p90stack[!grepl('slow', Framework) ]
i90stack = i90stack[!grepl('slow', Framework) ]


pdf(file.path(outdir, 'sim_scatter90.pdf'), width = 10, height = 5.5)
estscatter(rbind(estkw90himid,estkw90minhi), size=.4)
dev.off()

stop('scat')



point90r = sideside(p90stack, titl = '')
point90rzoom = sideside(p90stack[!(estimate %in% c('rev1', 'dm48') ), ], titl = '') 

ggsave(point90r, file = file.path(outdir, 'sim_rmse90.pdf'),
			 width = wid, height = hgt) 
ggsave(point90rzoom, file = file.path(outdir, 'sim_rmse90zoom.pdf'),
			 width = wid, height = hgt) 

point90bzoom = sideside(p90stack[estimate != 'dm48', ], metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = '')

ggsave(point90bzoom, file = file.path(outdir, 'sim_bias90zoom.pdf'),
			 width = wid, height = hgt) 

int90c = sideside(i90stack, metric = 'Coverage', titl = '', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100)
int90c = int90c + geom_hline(yintercept = c(85,95), lty=3)
int90w = sideside(i90stack, metric = 'Width', titl = '')
	
ggsave(int90c, file = file.path(outdir, 'sim_cover90.pdf'),
			 width = wid, height = hgt) 
ggsave(int90w, file = file.path(outdir, 'sim_width90.pdf'),
			 width = wid, height = hgt) 
	
cat(base::date(), '\n')
					
					



