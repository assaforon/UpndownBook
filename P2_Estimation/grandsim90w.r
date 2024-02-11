rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

logi90parm = fread(file.path(outdir, 'scenarios_logi90.csv') )
weib90parm = fread(file.path(outdir, 'scenarios_weib90.csv') )

### Comparative simulations for 90th percentile

#-------------------- Prep


#### Constants

nsim = 400 # This overrides 'nsim' in the scenario-generation scripts
M = 12
n = 60
ktarg90 = k2targ(6)
k90list = list(k=6, lowTarget=FALSE, fastStart=TRUE)
b90list = list(coin=1/9, lowTarget=FALSE, fastStart=TRUE)
lostart = 3
histart = 10

# Generate the scenarios on the dose grid
weib90F = weib90parm[1:nsim, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib90Fu = weib90parm[1:nsim, ][ ,as.vector(pweib3((1:M)-3, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")'
#  (we use the same thresholds throughout, to eliminate this variability source)
#  (but a separate one for each distribution, just in case we got a slightly odd draw)

thresh90w = matrix(runif(n*nsim), nrow=n)
thresh90l = matrix(runif(n*nsim), nrow=n)
					

	
#-------------------- simulation and estimation: Weibull	

########### k-row standard	

kw90lomid = dfsim(60, starting = lostart, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
estkw90lomid = estbatch(kw90lomid, truth=weib90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
		
kw90himid = dfsim(60, starting = histart, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
estkw90himid = estbatch(kw90himid, truth=weib90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)

kw90midmid = dfsim(60, starting = M/2, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
estkw90midmid = estbatch(kw90midmid, truth=weib90parm$t90[1:nsim], target=.9, bpt=ktarg90, 
				desfun=krow, desargs=k90list)

kw90minhi = dfsim(60, starting = 1, Fvals=weib90Fu, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
estkw90minhi = estbatch(kw90minhi, truth=weib90parm$t90[1:nsim]+3, target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
	
cat('k standard\n')	
### k-row variants

# Turning off fastStart
kw90himid_slow = dfsim(60, starting = histart, Fvals=weib90F, ensemble=nsim,
						desArgs=list(k=6, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90w)					
estkw90himid_slow = estbatch(kw90himid_slow, truth=weib90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
kw90lomid_slow = dfsim(60, starting = lostart, Fvals=weib90F, ensemble=nsim,
						desArgs=list(k=6, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90w)					
estkw90lomid_slow = estbatch(kw90lomid_slow, truth=weib90parm$t90[1:nsim], target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)

# Playing with k
kw90minhi5 = dfsim(60, starting = 1, Fvals=weib90Fu, ensemble=nsim,
						desArgs=list(k=5, lowTarget=FALSE, fastStart=TRUE), 
						thresholds = thresh90w)					
estkw90minhi5 = estbatch(kw90minhi5, truth=weib90parm$t90[1:nsim]+3, target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)
kw90minhi7 = dfsim(60, starting = 1, Fvals=weib90Fu, ensemble=nsim,
						desArgs=list(k=7, lowTarget=FALSE, fastStart=TRUE), 
						thresholds = thresh90w)					
estkw90minhi7 = estbatch(kw90minhi7, truth=weib90parm$t90[1:nsim]+3, target=.9, 
					bpt=ktarg90, desfun=krow, desargs=k90list)

cat('k variants\n')						
########### BCD standard	

bw90lomid = dfsim(60, starting = lostart, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90w)					
estbw90lomid = estbatch(bw90lomid, truth=weib90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)
		
bw90himid = dfsim(60, starting = histart, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90w)					
estbw90himid = estbatch(bw90himid, truth=weib90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)

bw90midmid = dfsim(60, starting = M/2, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90w)					
estbw90midmid = estbatch(bw90midmid, truth=weib90parm$t90[1:nsim], target=.9,  
				desfun=bcd, desargs=b90list)

bw90minhi = dfsim(60, starting = 1, Fvals=weib90Fu, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh90w)					
estbw90minhi = estbatch(bw90minhi, truth=weib90parm$t90[1:nsim]+3, target=.9, 
					 desfun=bcd, desargs=b90list)

cat('b standard\n')						
### BCD variants

# Turning off fastStart
bw90himid_slow = dfsim(60, starting = histart, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=list(coin=1/9, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90w)					
estbw90himid_slow = estbatch(bw90himid_slow, truth=weib90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)
bw90lomid_slow = dfsim(60, starting = lostart, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=list(coin=1/9, lowTarget=FALSE, fastStart=FALSE), 
						thresholds = thresh90w)					
estbw90lomid_slow = estbatch(bw90lomid_slow, truth=weib90parm$t90[1:nsim], target=.9, 
					 desfun=bcd, desargs=b90list)

save.image(file.path(outdir, 'grandsim90w.RData'))					
cat('b variants and done.\n')					
									





		

		
					





