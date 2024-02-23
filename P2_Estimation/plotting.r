cat(base::date(), '\n')
rm(list=ls())

#---------------------- Preamble
outdir = '../../output'

load(file.path(outdir, 'grandsim30w.RData'))	
load(file.path(outdir, 'grandsim30w_gud.RData'))	
e30w = c(ls(pat='est[bk]w'), ls(pat='estg[23]w') )
load(file.path(outdir, 'grandsim30l.RData'))	
e30l = c(ls(pat='est[bk]l3'), ls(pat='estg[23]l3') )

load(file.path(outdir, 'grandsim90w.RData'))	
e90w = ls(pat='est[bk]w9')
e90w = e90w[!grepl('slow',e90w) & !grepl('[57]$', e90w) ]

load(file.path(outdir, 'grandsim50w.RData'))	
load(file.path(outdir, 'grandsim50l.RData'))	
e50 = ls(pat='est[wl]5')


#### Constants

desnames = c('BCD', "K-row")
colors4 = c('grey65', 'grey40', 'grey80', 'black')
inames = c('Coverage', 'Width')

# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')


#---------------------- Munge

p90wstack = combo(e90w, finites = FALSE)
p90wstack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]
i90wstack = combo(e90w, atomfun = imetrix, outnames = inames)
i90wstack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]

p30wstack = combo(e30w, finites = FALSE)
p30wstack[ , Design := factor(substr(Framework, 1, 2), 
				labels = c('BCD', 'Group (2)', 'Group (3)', "K-row") ) ]
i30wstack = combo(e30w, atomfun = imetrix, outnames = inames)
# Designs; now with GUD
i30wstack[ , Design := factor(substr(Framework, 1, 2), 
				labels = c('BCD', 'Group (2)', 'Group (3)', "K-row") ) ]

p30lstack = combo(e30l, finites = FALSE)
p30lstack[ , Design := factor(substr(Framework, 1, 2), 
				labels = c('BCD', 'Group (2)', 'Group (3)', "K-row") ) ]
i30lstack = combo(e30l, atomfun = imetrix, outnames = inames)
i30lstack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]

p50stack = combo(e50, finites = FALSE)
p50stack[ , Design := factor(substr(Framework, 1, 1), 
			labels = c('Logistic', 'Weibull') )  ]
i50stack = combo(e50, atomfun = imetrix, outnames = inames)
i50stack[ , Design := factor(substr(Framework, 1, 1), 
			labels = c('Logistic', 'Weibull') )  ]


#---------------------- Plot

point90r = sideside(p90wstack, titl = 'RMSEs, 90th Percentile, Weibull Curves')
point90rzoom = sideside(p90wstack[estimate != 'dm48', ], titl = 'RMSEs, 90th Percentile, Weibull Curves, w/o Dixon-Mood') 
point90b = sideside(p90wstack, metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = 'Bias, 90th Percentile, Weibull Curves')
point90bzoom = sideside(p90wstack[estimate != 'dm48', ], metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = 'Bias, 90th Percentile, Weibull Curves, w/o Dixon-Mood')
point90qzoom = sideside(p90wstack[estimate != 'dm48', ], metric = 'QAE95', titl = 'QAE95, 90th Percentile, Weibull Curves, w/o Dixon-Mood') 

int90c = sideside(i90wstack, metric = 'Coverage', titl = 'Coverage, 90th Percentile, Weibull Curves', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100)
int90c = int90c + geom_hline(yintercept = c(85,95), lty=3)
int90w = sideside(i90wstack, metric = 'Width', titl = 'CI Width, 90th Percentile, Weibull Curves')

point30wr = sideside(p30wstack, titl = 'RMSEs, 30th Percentile, Weibull Curves', colkey=colors4)
point30wb = sideside(p30wstack, metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = 'Bias, 30th Percentile, Weibull Curves', colkey=colors4)

int30wc = sideside(i30wstack, metric = 'Coverage', titl = 'Coverage, 30th Percentile, Weibull Curves', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100, colkey=colors4)
int30wc = int30wc + geom_hline(yintercept = c(85,95), lty=3)
int30ww = sideside(i30wstack, metric = 'Width', titl = 'CI Width, 30th Percentile, Weibull Curves', colkey=colors4)

point30lr = sideside(p30lstack, titl = 'RMSEs, 30th Percentile, Logistic Curves')
point30lb = sideside(p30lstack, metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = 'Bias, 30th Percentile, Logistic Curves')

int30lc = sideside(i30lstack, metric = 'Coverage', titl = 'Coverage, 30th Percentile, Logistic Curves', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100)
int30lc = int30lc + geom_hline(yintercept = c(85,95), lty=3)
int30lw = sideside(i30lstack, metric = 'Width', titl = 'CI Width, 30th Percentile, Logistic Curves')

point50r = sideside(p50stack, titl = 'RMSEs, 50th Percentile')
point50r = point50r + labs(color='Family')
point50b = sideside(p50stack, metric = 'Bias', zoom=c(NA, NA), expansion = c(.01, .01), yref = 0, titl = 'Bias, 50th Percentile')
point50b = point50b + labs(color='Family')

int50c = sideside(i50stack, metric = 'Coverage', titl = 'Coverage, 50th Percentile', zoom=c(NA, NA), expansion = c(.01, .01), yref = 90, multip = 100)
int50c = int50c + geom_hline(yintercept = c(85,95), lty=3) + labs(color='Family')
int50w = sideside(i50stack, metric = 'Width', titl = 'CI Width, 50th Percentile')
int50w = int50w + labs(color='Family')


save.image(file.path(outdir, 'plotting.RData') )

cat(base::date(), '\n')










