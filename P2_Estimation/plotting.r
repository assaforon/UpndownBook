cat(base::date(), '\n')
rm(list=ls())

#---------------------- Preamble
outdir = '../../output'

load(file.path(outdir, 'grandsim30w.RData'))	
e30w = ls(pat='est[bk]w')

desnames = c('BCD', "K-row")
inames = c('Coverage', 'Width')

# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')


#---------------------- Munge

p30wstack = combo(e30w, finites = FALSE)
p30wstack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]
i30wstack = combo(e30w, atomfun = imetrix, outnames = inames)
i30wstack[ , Design := factor(substr(Framework, 1, 1), labels = desnames) ]


