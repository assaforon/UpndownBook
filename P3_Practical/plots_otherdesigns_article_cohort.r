cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(plyr)
library(ggplot2)
theme_set(theme_bw(19)) 
#library(patchwork)

load(file.path(outdir, 'cohort3othersim30w.RData'))	
all30w = ls(pat='w30[lmh][eio]')

#### Constants

wid1 = 11
hgt1 = 9

colors7 = c('orange', 'gold', 'turquoise4', 'turquoise1',  'turquoise3', 'dodgerblue',  'firebrick')
names6 = c('BOIN', 'CCD', paste('CRM', c('Low', 'Wide', 'High')), 'UDD (3,0,2)')
psize = 1.5
lsize = 1.3

source('summutils_other.r')

#--------------------------------- Main performance plots ----------------------------#

p30main = distill(all30w)
p30main[ , Design := mapvalues(des, sort(unique(des)), names6) ]
# Removing the "boring" CRM option
p30main = p30main[!grepl('Med', Design), ]

p30main2 = distill(all30w, yfun = function(x) sum(x$mtdest == x$mtd, na.rm = TRUE)/nrow(x) )
p30main2[ , Design := mapvalues(des, sort(unique(des)), names6) ]
p30main2 = p30main2[!grepl('Med', Design), ]

p30summ = p30main[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y) ), keyby = 'Design']
p30summ2 = p30main2[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y)), keyby = 'Design']


pm_30_0 <- ggplot(p30summ, aes(x, y, color = Design) ) + geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + 
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors7[-5]) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") 
pm_30 <- pm_30_0 + 
		coord_cartesian(ylim = c(60, 75), xlim = c(20,35), expand = 0) + geom_vline(xintercept = 30, lty = 2)

pm2_30 <- ggplot(p30summ2, aes(x, y, color = Design) ) + geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + 
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors7[-5]) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs that found the 'True MTD' (%)") #+ 
#		coord_cartesian(ylim = c(42.5, 57.5), xlim = c(20,35), expand = 0) + geom_vline(xintercept = 30, lty = 2)

#


#p1_30 <- ggplot(p30main, aes(100*x, 100*y, color = Design) ) + geom_point(size=psize) + scale_color_manual(values = colors7[-5]) +
#		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") + ylim(60, 80) + xlim(20,35) +
#		geom_vline(xintercept = 30, lty = 2)


ggsave(pm_30, file = file.path(outdir, 'othsim_main30_cohort.pdf'), width = wid1, height = hgt1)
ggsave(pm_30_0, file = file.path(outdir, 'othsim_main30_cohort_freelims.pdf'), width = wid1, height = hgt1)
ggsave(pm2_30, file = file.path(outdir, 'othsim_mtd30_cohort.pdf'), width = wid1, height = hgt1)
# ggsave(p1_30, file = file.path(outdir, 'othsim_main30.pdf'), width = wid1, height = hgt1)

#------------------------ "Number treated in window" histograms -----------------------#

e30w = ls(pat = 'rest[bcg][crow]')
e30w = e30w[grepl(30, e30w) & !grepl('4[_]05', e30w) & !grepl('boin', e30w) & !grepl('midmid', e30w) ]

e30combo = rbindlist(lapply(e30w, dextract), fill = TRUE)
e30combo[ , Design := mapvalues(des, sort(unique(des)), names6[-1] ) ]

e30combo[ , Setting := mapvalues(sett, c('minlo', 'minmid', 'minhi'), paste(c('Lower', 'Mid', 'Upper'), 'Target') ) ]

phist <- ggplot(e30combo, aes(ninterval/3)) + geom_histogram(fill='darkcyan', bins = 10) + facet_grid(Design ~ Setting) +
				labs(x = "Cohorts Treated in 'Acceptable Window'", y = "Number of Runs") + scale_x_continuous(breaks = seq(0,10,2))
				# + scale_x_continuous(limits=c(-1,31), expand = c(0,0))

ggsave(phist, file = file.path(outdir, 'othsim_hist30_cohort.pdf'), width = wid1, height = hgt1 * 1.4)
