cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

wid = 12
hgt = 7

load(file.path(outdir, 'grandsim50w.RData'))	
load(file.path(outdir, 'grandsim50l.RData'))	
e50 = ls(pat='est[wl]5')

#### Constants

desnames = c('BCD', "K-row")
inames = c('Coverage', 'Width')

gpoints = c("dm48", "dyna", "cir",  "ir",   "dual")

# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')


library(ggplot2)
theme_set(theme_bw(18))

p50stack = combo(e50)
p50stack[ , Design := factor(substr(Framework, 1, 1), 
			labels = c('Logistic', 'Weibull') )  ]
i50stack = combo(e50, atomfun = imetrix, outnames = inames, finites = FALSE)
i50stack[ , Design := factor(substr(Framework, 1, 1), 
			labels = c('Logistic', 'Weibull') )  ]

#------------------- Plots

point50r = sideside(p50stack)
point50r = point50r + labs(color='Family')
point50b = sideside(p50stack, metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0)
point50b = point50b + labs(color='Family')
point50q = sideside(p50stack, metric = 'MAE90') 
point50q = point50q + labs(color='Family')


int50c = sideside(i50stack, metric = 'Coverage', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100)
int50c = int50c + geom_hline(yintercept = c(85,95), lty=3) + labs(color='Family')
int50w = sideside(i50stack, metric = 'Width')
int50w = int50w + labs(color='Family')

stop(base::date(), 'Check plots!\n')



ggsave(point50r + labs(y = "RMSE (spacing units)", title = ''),
			file = file.path(outdir, 'sim_rmse50.pdf'), width=wid, height=hgt)

ggsave(point50b + labs(y = "Bias (spacing units)", title = ''),
			file = file.path(outdir, 'sim_bias50.pdf'), width=wid, height=hgt)

ggsave(int50c + labs(y = "Interval Coverage (%)", title = ''),
			file = file.path(outdir, 'sim_cover50.pdf'), width=wid, height=hgt)

ggsave(int50w + labs(y = "Average Interval Width (spacing units)", title = ''),
			file = file.path(outdir, 'sim_width50.pdf'), width=wid, height=hgt)
			
			
cat(base::date(), '\n')





