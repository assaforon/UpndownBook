rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

weib30parm = fread(file.path(outdir, 'scenarios_weib30.csv') )

### Comparative simulations for 30th percentile

#-------------------- Prep


#### Constants

set.seed(7928)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
M = 8
n = 60
ktarg30 = k2targ(2, lowTarget=TRUE)
gtarg302  = g2targ(3,0,2)
g2list = list(s=2, ll=0, ul=1)
g3list = list(s=3, ll=0, ul=2)
lostart = 1

# Generate the scenarios on the dose grid
weib30F = weib30parm[1:nsim, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib30Fu = weib30parm[1:nsim, ][ ,as.vector(pweib3((1:M)-2, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib30Fl = weib30parm[1:nsim, ][ ,as.vector(pweib3((1:M)+2, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh30w = matrix(runif(n*nsim), nrow=n)
					

	
#-------------------- simulation and estimation: Weibull	

########### 2 trial runs

g2w30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30w)					
estg2w30minmid = estbatch(g2w30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=groupUD, ccurvy = TRUE, desargs=g2list)
					
g3w30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30w)					
estg3w30minmid = estbatch(g3w30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
					 desfun=groupUD, ccurvy = TRUE, desargs=g3list, bpt = gtarg302)

	
# Good pausing point to see the output of one framework 
 
 #stop('ha')	

### cohort 2 runs

g2w30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30w)					
estg2w30midmid = estbatch(g2w30midmid, truth=weib30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=groupUD, ccurvy = TRUE, desargs=g2list)

g2w30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30w)					
estg2w30minhi = estbatch(g2w30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=groupUD, ccurvy = TRUE, desargs=g2list)

g2w30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30w)					
estg2w30minlo = estbatch(g2w30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=groupUD, ccurvy = TRUE, desargs=g2list)
	
cat('g2\n')	
					
########### cohort 3 runs

		
g3w30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30w)			
estg3w30midmid = estbatch(g3w30midmid, truth=weib30parm$t30[1:nsim], target=.3,  
				desfun=groupUD, ccurvy = TRUE, desargs=g3list, bpt = gtarg302)

g3w30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30w)			
estg3w30minhi = estbatch(g3w30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					 desfun=groupUD, ccurvy = TRUE, desargs=g3list, bpt = gtarg302)

g3w30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30w)			
estg3w30minlo = estbatch(g3w30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					 desfun=groupUD, ccurvy = TRUE, desargs=g3list, bpt = gtarg302)

				

save.image(file.path(outdir, 'grandsim30w_gud.RData'))					
cat('done.\n')					
									





		

		
					





