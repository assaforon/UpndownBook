rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

logi90parm = fread(file.path(outdir, 'scenarios_logi90.csv') )

### Comparative simulations for 90th percentile

#-------------------- Prep


#### Constants

set.seed(69942)
nsim = 100 # This overrides 'nsim' in the scenario-generation scripts
M = 12
n = 60
ktarg90 = k2targ(6)
k90list = list(k=6, lowTarget=FALSE, fastStart=TRUE)
b90list = list(coin=1/9, lowTarget=FALSE, fastStart=TRUE)
lostart = 3
histart = 10

# Generate the scenarios on the dose grid
logi90F = logi90parm[1:nsim, ][ ,as.vector(plogis(1:M, scale=scal, location=offs-shif)), by = 'row0'] %$% matrix(V1, nrow = M)

logi90Fu = logi90parm[1:nsim, ][ ,as.vector(plogis((1:M)-3, scale=scal, location=offs-shif)), by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")'
#  (we use the same thresholds throughout, to eliminate this variability source)
#  (but a separate one for each distribution, just in case we got a slightly odd draw)

thresh90l = matrix(runif(n*nsim), nrow=n)
					
	
#-------------------- simulation and estimation: Logistic	

########### k-row standard	

kl90lomid = dfsim(n, starting = lostart, Fvals=logi90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90l)					
estkl90lomid = estbatch(kl90lomid, truth=logi90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
		
kl90himid = dfsim(n, starting = histart, Fvals=logi90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90l)					
estkl90himid = estbatch(kl90himid, truth=logi90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)

kl90midmid = dfsim(n, starting = M/2, Fvals=logi90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90l)					
estkl90midmid = estbatch(kl90midmid, truth=logi90parm$t90[1:nsim], target=.9, bpt=ktarg90, 
				desfun=krow, desargs=k90list)

kl90minhi = dfsim(n, starting = 1, Fvals=logi90Fu, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90l)					
estkl90minhi = estbatch(kl90minhi, truth=logi90parm$t90[1:nsim]+3, target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
	
cat('k standard\n')	
### k-row variants

# Turning off fastStart
kl90himid_slow = dfsim(n, starting = histart, Fvals=logi90F, ensemble=nsim,
						desArgs=list(k=6, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90l)					
estkl90himid_slow = estbatch(kl90himid_slow, truth=logi90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
kl90lomid_slow = dfsim(n, starting = lostart, Fvals=logi90F, ensemble=nsim,
						desArgs=list(k=6, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90l)					
estkl90lomid_slow = estbatch(kl90lomid_slow, truth=logi90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)

# Playing with k
kl90minhi5 = dfsim(n, starting = 1, Fvals=logi90Fu, ensemble=nsim,
						desArgs=list(k=5, lowTarget=FALSE, fastStart=TRUE), 
						thresholds = thresh90l)					
estkl90minhi5 = estbatch(kl90minhi5, truth=logi90parm$t90[1:nsim]+3, target=.9, 
					bpt=k2targ(5), desfun=krow, desargs=k90list)
kl90minhi7 = dfsim(n, starting = 1, Fvals=logi90Fu, ensemble=nsim,
						desArgs=list(k=7, lowTarget=FALSE, fastStart=TRUE), 
						thresholds = thresh90l)					
estkl90minhi7 = estbatch(kl90minhi7, truth=logi90parm$t90[1:nsim]+3, target=.9, 
					bpt=k2targ(7), desfun=krow, desargs=k90list)

cat('k variants\n')						
########### BCD standard	

bl90lomid = dfsim(n, starting = lostart, Fvals=logi90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90l)					
estbl90lomid = estbatch(bl90lomid, truth=logi90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)
		
bl90himid = dfsim(n, starting = histart, Fvals=logi90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90l)					
estbl90himid = estbatch(bl90himid, truth=logi90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)

bl90midmid = dfsim(n, starting = M/2, Fvals=logi90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90l)					
estbl90midmid = estbatch(bl90midmid, truth=logi90parm$t90[1:nsim], target=.9,  
				desfun=bcd, desargs=b90list)

bl90minhi = dfsim(n, starting = 1, Fvals=logi90Fu, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90l)					
estbl90minhi = estbatch(bl90minhi, truth=logi90parm$t90[1:nsim]+3, target=.9, 
					 desfun=bcd, desargs=b90list)

cat('b standard\n')						
### BCD variants

# Turning off fastStart
bl90himid_slow = dfsim(n, starting = histart, Fvals=logi90F, ensemble=nsim,
						design=bcd, desArgs=list(coin=1/9, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90l)					
estbl90himid_slow = estbatch(bl90himid_slow, truth=logi90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)
bl90lomid_slow = dfsim(n, starting = lostart, Fvals=logi90F, ensemble=nsim,
						design=bcd, desArgs=list(coin=1/9, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90l)					
estbl90lomid_slow = estbatch(bl90lomid_slow, truth=logi90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)

save.image(file.path(outdir, 'grandsim90l.RData'))					
cat('b variants and done.\n')					
									





		

		
					





