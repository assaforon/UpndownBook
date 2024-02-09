rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

logi90parm = fread(file.path(outdir, 'scenarios_logi90.csv') )
weib90parm = fread(file.path(outdir, 'scenarios_weib90.csv') )

#### Constants

nsim = 500
M = 12
n = 60
ktarg90 = k2targ(6)
k90list = list(k=6, lowTarget=FALSE, fastStart=TRUE)

weib90F = weib90parm[1:nsim, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
w90lomid = dfsim(60, starting = 4, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list)					

w90lomid_est = estbatch(w90lomid, truth=weib90parm$t90[1:nsim], target=.9, bpt=ktarg90, 
				desfun=krow, desargs=k90list)
				# , addLiao = TRUE) 

		

		
					





