cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(plyr)
library(ggplot2)
theme_set(theme_bw(19)) 
#library(patchwork)

load(file.path(outdir, 'othersim30w.RData'))	
all30w = ls(pat='w30[lmh][eio]')
load(file.path(outdir, 'othersim90w.RData'))	
all90w = ls(pat='w90[lmh][eio]')

#### Constants

wid1 = 11
hgt1 = 9

colors7 = c('orange', 'gold', 'turquoise4', 'turquoise1',  'turquoise3', 'dodgerblue',  'firebrick')
names7 = c('BOIN', 'CCD', paste('CRM', c('Low', 'Med', 'Wide', 'High')), 'UDD (k=2)')
colors6 = c('gold', 'turquoise4', 'turquoise1',  'turquoise3', 'dodgerblue', 'firebrick')
names6 = c('CCD', paste('CRM', c('High', 'Wide', 'Low', 'Med')), 'UDD (k=6)')
psize = 1.5
lsize = 1.3

source('summutils_other.r')

#--------------------------------- Main performance plots ----------------------------#

p30main = distill(all30w)
p30main[ , Design := mapvalues(des, sort(unique(des)), names7) ]
# Removing the "boring" CRM option
p30main = p30main[!grepl('Med', Design), ]

p30main2 = distill(all30w, yfun = function(x) sum(x$mtdest == x$mtd, na.rm = TRUE)/nrow(x) )
p30main2[ , Design := mapvalues(des, sort(unique(des)), names7) ]
p30main2 = p30main2[!grepl('Med', Design), ]

p30summ = p30main[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y) ), keyby = 'Design']
p30summ2 = p30main2[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y)), keyby = 'Design']


pm_30 <- ggplot(p30summ, aes(x, y, color = Design) ) + geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + 
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors7[-5]) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") + 
		coord_cartesian(ylim = c(60, 75), xlim = c(20,35), expand = 0) + geom_vline(xintercept = 30, lty = 2)

pm2_30 <- ggplot(p30summ2, aes(x, y, color = Design) ) + geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + 
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors7[-5]) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs that found the 'True MTD' (%)") + 
		coord_cartesian(ylim = c(42.5, 57.5), xlim = c(20,35), expand = 0) + geom_vline(xintercept = 30, lty = 2)

#


#p1_30 <- ggplot(p30main, aes(100*x, 100*y, color = Design) ) + geom_point(size=psize) + scale_color_manual(values = colors7[-5]) +
#		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") + ylim(60, 80) + xlim(20,35) +
#		geom_vline(xintercept = 30, lty = 2)

p90main = distill(all90w)

p90main = p90main[des != 'ccd', ]
p90main[ , Design := mapvalues(des, sort(unique(des)), names6) ]
# Removing the "boring" CRM option
p90main = p90main[!grepl('Med', Design), ]
p90summ = p90main[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y) ), keyby = 'Design']

pm_90 <- ggplot(p90summ, aes(x, y, color = Design) ) + geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + 
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors6[-4]) +
		labs(x = "Ensemble-Mean Efficacy Rate (%)", y = "Runs with 'Best Dose' Estimate in 'Desirable Window' (%)")   + 
		coord_cartesian(ylim = c(58, 82), xlim = c(68,95), expand = 0) + geom_vline(xintercept = 90, lty = 2)



# p1_90 <- ggplot(p90main, aes(100*x, 100*y, color = Design) ) + geom_point(size=psize) + scale_color_manual(values = colors6[-4]) +
		# labs(x = "Ensemble-Mean Efficacy Rate (%)", y = "Runs with 'Best Dose' Estimate in 'Desirable Window' (%)") + ylim(50, 89) + xlim(65,95) +
		# geom_vline(xintercept = 90, lty = 2)

#### Saving

fwrite(p30main, file = file.path(outdir, 'sim30summaries.csv') )
fwrite(p30main2, file = file.path(outdir, 'sim30summaries_MTD.csv') )
fwrite(p90main, file = file.path(outdir, 'sim90summaries.csv') )

ggsave(pm_30, file = file.path(outdir, 'othsim_main30_RnR.pdf'), width = wid1, height = hgt1)
ggsave(pm2_30, file = file.path(outdir, 'othsim_mtd30_RnR.pdf'), width = wid1, height = hgt1)
ggsave(pm_90, file = file.path(outdir, 'othsim_main90_RnR.pdf'), width = wid1, height = hgt1)

# ggsave(p1_30, file = file.path(outdir, 'othsim_main30.pdf'), width = wid1, height = hgt1)
# ggsave(p1_90, file = file.path(outdir, 'othsim_main90.pdf'), width = wid1, height = hgt1)

#------------------------ "Number treated in window" histograms -----------------------#

#stop('Revision')

e30w = ls(pat = 'rest[bck][crow]')
e30w = e30w[grepl(30, e30w) & !grepl('4[_]05', e30w) & !grepl('boin', e30w) & !grepl('midmid', e30w) ]

e30combo = rbindlist(lapply(e30w, dextract), fill = TRUE)
e30combo[ , Design := mapvalues(des, sort(unique(des)), names7[c(2,3,5:7)] ) ]

e30combo[ , Setting := mapvalues(sett, c('minlo', 'minmid', 'minhi'), paste(c('Lower', 'Mid', 'Upper'), 'Target') ) ]

phist <- ggplot(e30combo, aes(ninterval)) + geom_histogram(fill='steelblue2', bins = 31) + facet_grid(Design ~ Setting) +
				labs(x = "Patients Treated in 'Acceptable Window'", y = "Number of Runs") 
				# + scale_x_continuous(limits=c(-1,31), expand = c(0,0))

ggsave(phist, file = file.path(outdir, 'othsim_hist30.pdf'), width = wid1, height = hgt1 * 1.4)
