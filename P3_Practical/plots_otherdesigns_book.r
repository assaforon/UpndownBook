cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

library(plyr)
library(ggplot2)
theme_set(theme_bw(17)) 
#library(patchwork)

load(file.path(outdir, 'othersim30w.RData'))	
all30w = ls(pat='w30[lmh][eio]')
all30w = all30w[!grepl('4[_]05',all30w)]
load(file.path(outdir, 'othersim90w.RData'))	
all90w = ls(pat='w90[lmh][eio]')
all90w = all90w[!grepl('7[_]025',all90w)]

#### Constants

wid1 = 11
hgt1 = 9
wid2 = 12
hgt2 = 7


colors6 = c(rep('grey75', 2), rep('grey50', 3), 'black')
names6 = c('BOIN', 'CCD', paste('CRM', c('Low', 'Wide', 'High')), 'UDD (k=2)')
shapes6 = c(18, 12, 25, 17, 15, 19)
#colors6 = c('gold', 'turquoise4', 'turquoise1',  'turquoise3', 'dodgerblue', 'firebrick')
#names6 = c('CCD', paste('CRM', c('High', 'Wide', 'Low', 'Med')), 'UDD (k=6)')
psize = 1.5
lsize = 1.3


source('summutils_other.r')

#--------------------------------- Main performance plots ----------------------------#

p30main = distill(all30w)
p30main[ , Design := mapvalues(des, sort(unique(des)), names6) ]
# Removing the "boring" CRM option
p30main = p30main[!grepl('Med', Design), ]

# stop('distilled 30')

p30summ = p30main[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y) ), keyby = 'Design']


pm_30 <- ggplot(p30summ, aes(x, y, color = Design, shape = Design, fill = Design) ) + 
		geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + scale_shape_manual(values = shapes6) +
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors6) +
		labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") + 
		coord_cartesian(ylim = c(60, 75), xlim = c(20,35), expand = 0) + geom_vline(xintercept = 30, lty = 2) +
		scale_fill_manual(values = colors6) 


# p1_30 <- ggplot(p30main, aes(100*x, 100*y, color = Design, shape = Design) ) + geom_point(size=psize) + scale_color_manual(values = colors6) +
		# labs(x = "Ensemble-Mean DLT Rate (%)", y = "Runs with MTD Estimate in 'Acceptable Window' (%)") + ylim(60, 76) + xlim(21,34) +
		# geom_vline(xintercept = 30, lty = 2) + scale_shape_manual(values = shapes6)
		
p90main = distill(all90w)
p90main[ , Design := mapvalues(des, sort(unique(des)), gsub(2, 6, names6[c(2, 2, 5:3, 6) ]) ) ]
p90main[des=='ccd', Design := 'CCD (0.09)' ]
p90main[des=='ccd 2', Design := 'CCD (0.05)' ]

# Removing the "boring" CRM option
p90main = p90main[!grepl('Med', Design), ]

p90summ = p90main[ , list(x = 100*mean(x), xmin = 100*min(x), xmax = 100*max(x), y = 100*mean(y), ymin = 100*min(y), ymax = 100*max(y) ), keyby = 'Design']


pm_90 <- ggplot(p90summ, aes(x, y, color = Design, shape = Design, fill = Design) ) + 
		geom_pointrange(aes(ymin=ymin, ymax=ymax), size=psize, lwd = lsize) + scale_shape_manual(values = c(1, shapes6[-1]) ) +
		geom_linerange(aes(xmin=xmin, xmax=xmax), lwd=lsize ) + scale_color_manual(values = colors6) +
		labs(x = "Ensemble-Mean Efficacy Rate (%)", y = "Runs with 'Best Dose' Estimate in 'Desirable Window' (%)") + 
		coord_cartesian(ylim = c(62, 82), xlim = c(68,95), expand = 0) + geom_vline(xintercept = 90, lty = 2) +
		scale_fill_manual(values = colors6) 

# p1_90 <- ggplot(p90main, aes(100*x, 100*y, color = Design, shape = Design) ) + geom_point(size=psize) + scale_color_manual(values = c(colors6[1], colors6[-1]) ) +
		# labs(x = "Ensemble-Mean Efficacy Rate (%)", y = "Runs with 'Best Dose' Estimate in 'Desirable Window' (%)") + ylim(50, 89) + xlim(65,95) +
		# geom_vline(xintercept = 90, lty = 2) + scale_shape_manual(values = shapes6[c(6, 2:6) ])

ggsave(pm_30, file = file.path(outdir, 'othsim_main30_book.pdf'), width = wid1, height = hgt1)
ggsave(pm_90, file = file.path(outdir, 'othsim_main90_book.pdf'), width = wid1, height = hgt1)

# stop('enuff')

#------------------------ "Number treated in window" histograms -----------------------#

e30w0 = ls(pat = 'rest[bck][crow]')
e30w = e30w0[grepl(30, e30w0) & !grepl('4[_]05', e30w0) & !grepl('boin', e30w0) & !grepl('midmid', e30w0) ]

e30combo = rbindlist(lapply(e30w, dextract), fill = TRUE)
e30combo[ , Design := p30main$Design[match(des,p30main$des)]  ]

e30combo[ , Setting := mapvalues(sett, c('minlo', 'minmid', 'minhi'), paste(c('Lower', 'Mid', 'Upper'), 'Target') ) ]

phist <- ggplot(e30combo, aes(ninterval)) + geom_histogram(bins = 31) + facet_grid(Design ~ Setting) +
				labs(x = "Patients Treated in 'Acceptable Window'", y = "Number of Runs") 
				# + scale_x_continuous(limits=c(-1,31), expand = c(0,0))

ggsave(phist, file = file.path(outdir, 'othsim_hist30_book.pdf'), width = wid1, height = hgt1 * 1.4)

cat(base::date(), 'Article grayscaled\n')

#---------------------------- Continuous metrics

jwid = 0.1
theme_set(theme_bw(18)) 
# source('../P2_Estimation/sumsims.r')

e30w2 = e30w0[grepl(30, e30w0) & !grepl('4[_]05', e30w0)] 
e30pcombo = rbindlist(lapply(e30w2, dextract), fill = TRUE)

pmetrix30 = e30pcombo[ , list(RMSE = rmse(pointest, true, winsor = TRUE), MAE90 = mae(pointest, true),
								Bias = bias(pointest, true) ), keyby = .(des, sett)] 
pmetrix30[ , Design := mapvalues(des, sort(unique(des)), names6) ]

e90w2 = e30w0[grepl(90, e30w0) & !grepl('7[_]025', e30w0)] 
e90pcombo = rbindlist(lapply(e90w2, dextract), fill = TRUE)

pmetrix90 = e90pcombo[ , list(RMSE = rmse(pointest, true, winsor = TRUE), MAE90 = mae(pointest, true),
								Bias = bias(pointest, true)), keyby = .(des, sett)] 
pmetrix90[ , Design := mapvalues(des, sort(unique(des)), gsub(2, 6, names6[c(2, 2, 5:3, 6) ]) ) ]
pmetrix90[des=='ccd', Design := 'CCD (0.09)' ]
pmetrix90[des=='ccd 2', Design := 'CCD (0.05)' ]

cont30m = ggplot(pmetrix30, aes(Design, MAE90)) + geom_jitter(size = psize-1, width = jwid) + ylim(0,NA) + xlab('')
cont30r = ggplot(pmetrix30, aes(Design, RMSE)) + geom_jitter(size = psize-1, width = jwid) + ylim(0,NA) + xlab('')
cont30b = ggplot(pmetrix30, aes(Design, Bias)) + geom_jitter(size = psize-1, width = jwid) + geom_hline(yintercept = 0) + xlab('')

ggsave(cont30m, file = file.path(outdir, 'othsim30_mae.pdf'), width = wid2, height = hgt2)
ggsave(cont30r, file = file.path(outdir, 'othsim30_rmse.pdf'), width = wid2, height = hgt2)
ggsave(cont30b, file = file.path(outdir, 'othsim30_bias.pdf'), width = wid2, height = hgt2)

cont90m = ggplot(pmetrix90, aes(Design, MAE90)) + geom_jitter(size = psize-1, width = jwid) + ylim(0,NA) + xlab('')
cont90r = ggplot(pmetrix90, aes(Design, RMSE)) + geom_jitter(size = psize-1, width = jwid) + ylim(0,NA) + xlab('')
cont90b = ggplot(pmetrix90, aes(Design, Bias)) + geom_jitter(size = psize-1, width = jwid) + geom_hline(yintercept = 0) + xlab('')

ggsave(cont90m, file = file.path(outdir, 'othsim90_mae.pdf'), width = wid2, height = hgt2)
ggsave(cont90r, file = file.path(outdir, 'othsim90_rmse.pdf'), width = wid2, height = hgt2)
ggsave(cont90b, file = file.path(outdir, 'othsim90_bias.pdf'), width = wid2, height = hgt2)

cat(base::date(), 'Done.\n')

