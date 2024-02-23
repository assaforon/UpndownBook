rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

logi30parm = fread(file.path(outdir, 'scenarios_logi30.csv') )

### Comparative simulations for 30th percentile

#-------------------- Prep


#### Constants

set.seed(792333)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
M = 8
n = 60
ktarg30 = k2targ(2, lowTarget=TRUE)
gtarg302  = g2targ(3,0,2)
g2list = list(s=2, ll=0, ul=1)
g3list = list(s=3, ll=0, ul=2)
lostart = 1

# Generate the scenarios on the dose grid
logi30F = logi30parm[1:nsim, ][ ,as.vector(plogis(1:M, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
logi30Fu = logi30parm[1:nsim, ][ ,as.vector(plogis((1:M)-2, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
logi30Fl = logi30parm[1:nsim, ][ ,as.vector(plogis((1:M)+2, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh30l = matrix(runif(n*nsim), nrow=n)
					

	
#-------------------- simulation and estimation: logiull	

########### 2 trial runs

g2l30minmid = dfsim(n, starting = lostart, Fvals=logi30F, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30l)					
estg2l30minmid = estbatch(g2l30minmid, truth=logi30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=groupUD, desargs=g2list)
					
g3l30minmid = dfsim(n, starting = lostart, Fvals=logi30F, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30l)					
estg3l30minmid = estbatch(g3l30minmid, truth=logi30parm$t30[1:nsim], target=.3, 
					 desfun=groupUD, desargs=g3list, bpt = gtarg302)

	
# Good pausing point to see the output of one framework 
 
 #stop('ha')	

### cohort 2 runs

g2l30midmid = dfsim(n, starting = M/2, Fvals=logi30F, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30l)					
estg2l30midmid = estbatch(g2l30midmid, truth=logi30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=groupUD, desargs=g2list)

g2l30minhi = dfsim(n, starting = 1, Fvals=logi30Fu, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30l)					
estg2l30minhi = estbatch(g2l30minhi, truth=logi30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=groupUD, desargs=g2list)

g2l30minlo = dfsim(n, starting = 1, Fvals=logi30Fl, ensemble=nsim, design=groupUD,
						desArgs=g2list, thresholds = thresh30l)					
estg2l30minlo = estbatch(g2l30minlo, truth=logi30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=groupUD, desargs=g2list)
	
cat('g2\n')	
					
########### cohort 3 runs

		
g3l30midmid = dfsim(n, starting = M/2, Fvals=logi30F, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30l)			
estg3l30midmid = estbatch(g3l30midmid, truth=logi30parm$t30[1:nsim], target=.3,  
				desfun=groupUD, desargs=g3list, bpt = gtarg302)

g3l30minhi = dfsim(n, starting = 1, Fvals=logi30Fu, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30l)			
estg3l30minhi = estbatch(g3l30minhi, truth=logi30parm$t30[1:nsim]+2, target=.3, 
					 desfun=groupUD, desargs=g3list, bpt = gtarg302)

g3l30minlo = dfsim(n, starting = 1, Fvals=logi30Fl, ensemble=nsim,
						design=groupUD, desArgs=g3list, thresholds = thresh30l)			
estg3l30minlo = estbatch(g3l30minlo, truth=logi30parm$t30[1:nsim]-2, target=.3, 
					 desfun=groupUD, desargs=g3list, bpt = gtarg302)

				

save.image(file.path(outdir, 'grandsim30l_gud.RData'))					
cat('done.\n')					
									





		

		
					





