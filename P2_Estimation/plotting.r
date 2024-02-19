cat(base::date(), '\n')
rm(list=ls())

#---------------------- Preamble
source('simulation_header.r')

load(file.path(outdir, 'grandsim30w.RData'))	
e30w = ls(pat='est[bk]w')


# sourcing this code last, in case datasets include older copies of functions
source('sumsims.r')

p30wstack = combo(e30w, finites = FALSE)
i30wstack = combo(e30w, atomfun = imetrix, outnames = c('Coverage', 'Width'))


