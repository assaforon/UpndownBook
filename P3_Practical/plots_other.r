cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(plyr)
library(ggplot2)
theme_set(theme_bw(20)) 
#library(patchwork)

load(file.path(outdir, 'othersim30w.RData'))	
all30w = ls(pat='w30[lmh][eio]')
load(file.path(outdir, 'othersim90w.RData'))	
all90w = ls(pat='w90[lmh][eio]')

#### Constants

wid1 = 9
hgt1 = 9

colors7 = c('orange', 'gold', 'turquoise4', 'turquoise1',  'turquoise3', 'blue', 'firebrick')
names7 = c('BOIN', 'CCD', paste('CRM', c('Low', 'Med', 'Wide', 'High')), 'K-Row')
colors6 = c('gold', 'turquoise4', 'turquoise1',  'turquoise3', 'blue', 'firebrick')
names6 = c('CCD', paste('CRM', c('High', 'Wide', 'Low', 'Med')), 'K-Row')


source('summutils_other.r')

#--------------------------------- Main performance plots ----------------------------#

p30main = distill(all30w)
p30main[ , Design := mapvalues(des, sort(unique(des)), names7) ]

p1_30 <- ggplot(p30main, aes(100*x, 100*y, color = Design) ) + geom_point(size=4) + scale_color_manual(values = colors7) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Good Window' (%)") + ylim(50, 89) + xlim(20,35) +
		geom_vline(xintercept = 30, lty = 2)

p90main = distill(all90w)
p90main[ , Design := mapvalues(des, sort(unique(des)), names6) ]

p1_90 <- ggplot(p90main, aes(100*x, 100*y, color = Design) ) + geom_point(size=4) + scale_color_manual(values = colors6) +
		labs(x = "Ensemble-Mean Efficacy Rate (%)", y = "Runs with 'Best Dose' Estimate in 'Good Window' (%)") + ylim(50, 89) + xlim(65,95) +
		geom_vline(xintercept = 90, lty = 2)

ggsave(p1_30, file = file.path(outdir, 'othsim_main30.pdf'), width = wid1, height = hgt1)
ggsave(p1_90, file = file.path(outdir, 'othsim_main90.pdf'), width = wid1, height = hgt1)

#------------------------ "Number treated in window" histograms -----------------------#

e30w = ls(pat = 'rest[bck][crow]')
e30w = e30w[grepl(30, e30w) & !grepl('4[_]05', e30w) & !grepl('boin', e30w) & !grepl('midmid', e30w) ]

e30combo = rbindlist(lapply(e30w, dextract), fill = TRUE)
e30combo[ , Design := mapvalues(des, sort(unique(des)), names7[c(2,3,5:7)] ) ]

ggplot(e30combo, aes(ninterval)) + geom_histogram() + facet_grid(

