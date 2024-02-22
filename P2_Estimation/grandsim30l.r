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
k30list = list(k=2, lowTarget=TRUE, fastStart=FALSE)
b30list = list(coin=3/7, lowTarget=TRUE, fastStart=FALSE)
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
					
#-------------------- simulation and estimation: logistic

########### k-row standard	

kl30minmid = dfsim(n, starting = lostart, Fvals=logi30F, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30l)					
estkl30minmid = estbatch(kl30minmid, truth=logi30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
	
# Good pausing point to see the output of one framework 
 
# stop('ha')	

kl30midmid = dfsim(n, starting = M/2, Fvals=logi30F, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30l)					
estkl30midmid = estbatch(kl30midmid, truth=logi30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=krow, desargs=k30list)

kl30minhi = dfsim(n, starting = 1, Fvals=logi30Fu, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30l)					
estkl30minhi = estbatch(kl30minhi, truth=logi30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
# estkl30minhi_2 = estbatch(kl30minhi, truth=logi30parm$t30[1:nsim]+2, target=.3, 
	#				bpt=ktarg30, desfun=krow, desargs=k30list, randboot = FALSE)

kl30minlo = dfsim(n, starting = 1, Fvals=logi30Fl, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30l)					
estkl30minlo = estbatch(kl30minlo, truth=logi30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
	
cat('k standard\n')	
					
########### BCD standard	

bl30minmid = dfsim(n, starting = lostart, Fvals=logi30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30l)					
estbl30minmid = estbatch(bl30minmid, truth=logi30parm$t30[1:nsim], target=.3, 
					 desfun=bcd, desargs=b30list)
		
bl30midmid = dfsim(n, starting = M/2, Fvals=logi30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30l)			
estbl30midmid = estbatch(bl30midmid, truth=logi30parm$t30[1:nsim], target=.3,  
				desfun=bcd, desargs=b30list)

bl30minhi = dfsim(n, starting = 1, Fvals=logi30Fu, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30l)			
estbl30minhi = estbatch(bl30minhi, truth=logi30parm$t30[1:nsim]+2, target=.3, 
					 desfun=bcd, desargs=b30list)

bl30minlo = dfsim(n, starting = 1, Fvals=logi30Fl, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30l)			
estbl30minlo = estbatch(bl30minlo, truth=logi30parm$t30[1:nsim]-2, target=.3, 
					 desfun=bcd, desargs=b30list)

cat('b standard\n')						

save.image(file.path(outdir, 'grandsim30l.RData'))					
cat('done.\n')					
									





		

		
					





