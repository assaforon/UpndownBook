rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

logi50parm = fread(file.path(outdir, 'scenarios_logi50.csv') )

### Comparative simulations for 50th percentile - Logistic curves (symmetric around target)

#-------------------- Prep


#### Constants

set.seed(79280)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
M = 12
n = 40
lostart = 3
histart = 10

# Generate the scenarios on the dose grid
logi50F = logi50parm[1:nsim, ][ ,as.vector(plogis(1:M, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
logi50Fu = logi50parm[1:nsim, ][ ,as.vector(plogis((1:M)-3, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
logi50Fl = logi50parm[1:nsim, ][ ,as.vector(plogis((1:M)+3, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh50l = matrix(runif(n*nsim), nrow=n)
					

	
#-------------------- simulation and estimation	

l50lomid = dfsim(n, starting = lostart, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l)					
estl50lomid = estbatch(l50lomid, truth=logi50parm$t50[1:nsim], target=.5)
	
# Good pausing point to see the output of one framework 
# stop('ha')	

l50himid = dfsim(n, starting = histart, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l)					
estl50himid = estbatch(l50himid, truth=logi50parm$t50[1:nsim], target=.5)

l50midmid = dfsim(n, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l)					
estl50midmid = estbatch(l50midmid, truth=logi50parm$t50[1:nsim], target=.5)

l50minhi = dfsim(n, starting = 1, Fvals=logi50Fu, ensemble=nsim,
						thresholds = thresh50l)					
estl50minhi = estbatch(l50minhi, truth=logi50parm$t50[1:nsim]+3, target=.5)
	
l50maxlo = dfsim(n, starting = M, Fvals=logi50Fl, ensemble=nsim,
						thresholds = thresh50l)					
estl50maxlo = estbatch(l50maxlo, truth=logi50parm$t50[1:nsim]-3, target=.5)
				

save.image(file.path(outdir, 'grandsim50l.RData'))					
cat('done.\n')					
									





		

		
					





