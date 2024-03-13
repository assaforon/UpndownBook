cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(ggplot2)
theme_set(theme_bw(18)) 
library(patchwork)

load(file.path(outdir, 'halfsim30w.RData'))	
e30w = ls(pat='[h]*est[bgk][23]*w3')

#### Constants

wid = 12
hgt = 7

des2 = c('BCD', "K-row")
des4 = c('BCD', 'Group (2)', 'Group (3)', "K-row")
colors2 = c('grey65', 'black')
colors4 = c('grey65', 'grey40', 'grey80', 'black')
inames = c('Coverage', 'Width')

gpoints = c("dm48", "dyna", "cir",  "ir",   "dual")

# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')

p30stack = combo(e30w)
p30stack[ , Design := factor(substr(gsub('^h','',Framework), 1, 2), labels = des4) ]
i30stack = combo(e30w, atomfun = imetrix, outnames = inames, finites = FALSE)
i30stack[ , Design := factor(substr(gsub('^h','',Framework), 1, 2), labels = des4) ]


point30n1r = sideside(p30stack[!grepl('Group',Design) & grepl('^h',Framework), ])
ggsave(point30n1r, file = file.path(outdir, 'sim_rmse30n1.pdf'),
			 width = wid, height = hgt) 

point30n1rg = sideside(p30stack[estimate %in% gpoints & grepl('^h',Framework), ],
					colkey = colors4)
ggsave(point30n1rg, file = file.path(outdir, 'sim_rmse30n1g.pdf'),
			 width = wid, height = hgt) 

#sideside(p30stack[estimate %in% gpoints & !grepl('^h',Framework), ],metric='MAE90', colkey = colors4)

#sideside(p30stack[estimate %in% gpoints & !grepl('^h',Framework), ],metric='Bias', colkey = colors4,zoom=c(NA, NA), expansion = c(.01, .01), yref = 0)



