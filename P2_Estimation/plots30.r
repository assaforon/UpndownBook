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


#-------- Plots: point, n=30

point30n1r = sideside(p30stack[grepl('^h',Framework), ], colkey = colors4)
ggsave(point30n1r, file = file.path(outdir, 'sim_rmse30n1.pdf'),
			 width = wid, height = hgt) 


point30n1b = sideside(p30stack[grepl('^h',Framework),],
	metric = 'Bias', colkey = colors4, zoom=c(NA, NA), expansion = c(.01, .01), yref = 0)
ggsave(point30n1b, file = file.path(outdir, 'sim_bias30n1.pdf'),
			 width = wid, height = hgt) 

#-------- Plots: point, n=60

point30n2r = sideside(p30stack[!grepl('^h',Framework), ], colkey = colors4)
ggsave(point30n2r, file = file.path(outdir, 'sim_rmse30n2.pdf'),
			 width = wid, height = hgt) 

point30n2b = sideside(p30stack[!grepl('^h',Framework),],
	metric = 'Bias', colkey = colors4, zoom=c(NA, NA), expansion = c(.01, .01), yref = 0)
ggsave(point30n2b, file = file.path(outdir, 'sim_bias30n2.pdf'),
			 width = wid, height = hgt) 


#-------- Plots: interval

int30n1c = sideside(i30stack[grepl('^h',Framework), ], metric = 'Coverage', titl = '', 
		zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100, colkey = colors4)
ggsave(int30n1c + geom_hline(yintercept = c(85,95), lty=3) + labs(y = 'Interval Coverage (%)') , file = file.path(outdir, 'sim_cover30n1.pdf'),
			 width = wid, height = hgt) 

int30n2c = sideside(i30stack[!grepl('^h',Framework), ], metric = 'Coverage', titl = '', 
	zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100, colkey = colors4)
ggsave(int30n2c + geom_hline(yintercept = c(85,95), lty=3) + labs(y = 'Interval Coverage (%)') , file = file.path(outdir, 'sim_cover30n2.pdf'),
			 width = wid, height = hgt) 


int30n1w = sideside(i30stack[grepl('^h',Framework), ], metric = 'Width', titl = '', colkey = colors4)
ggsave(int30n1w, file = file.path(outdir, 'sim_width30n1.pdf'),
			 width = wid, height = hgt) 

int30n2w = sideside(i30stack[!grepl('^h',Framework), ], metric = 'Width', titl = '', colkey = colors4)
ggsave(int30n2w, file = file.path(outdir, 'sim_width30n2.pdf'),
			 width = wid, height = hgt) 



#sideside(p30stack[estimate %in% gpoints & !grepl('^h',Framework), ],metric='MAE90', colkey = colors4)

#sideside(p30stack[estimate %in% gpoints & !grepl('^h',Framework), ],metric='Bias', colkey = colors4,zoom=c(NA, NA), expansion = c(.01, .01), yref = 0)



