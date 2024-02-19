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
n = 50
ktarg30 = k2targ(2, lowTarget=TRUE)
k30list = list(k=2, lowTarget=TRUE, fastStart=FALSE)
b30list = list(coin=3/7, lowTarget=TRUE, fastStart=FALSE)
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

########### k-row standard	

kw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30minmid = estbatch(kw30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
	
# Good pausing point to see the output of one framework 
 
# stop('ha')	

kw30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30midmid = estbatch(kw30midmid, truth=weib30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=krow, desargs=k30list)

kw30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30minhi = estbatch(kw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
# estkw30minhi_2 = estbatch(kw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
	#				bpt=ktarg30, desfun=krow, desargs=k30list, randboot = FALSE)

kw30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30minlo = estbatch(kw30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=krow, desargs=k30list)
	
cat('k standard\n')	
					
########### BCD standard	

bw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)					
estbw30minmid = estbatch(bw30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
					 desfun=bcd, desargs=b30list)
		
bw30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30midmid = estbatch(bw30midmid, truth=weib30parm$t30[1:nsim], target=.3,  
				desfun=bcd, desargs=b30list)

bw30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30minhi = estbatch(bw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					 desfun=bcd, desargs=b30list)

bw30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30minlo = estbatch(bw30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					 desfun=bcd, desargs=b30list)

cat('b standard\n')						

save.image(file.path(outdir, 'grandsim30w.RData'))					
cat('b variants and done.\n')					
									





		

		
					





