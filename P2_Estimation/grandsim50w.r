rm(list=ls())
library(magrittr)

#---------------------- Preamble

source('simulation_header.r')
cat(base::date(), '\n')

weib50parm = fread(file.path(outdir, 'scenarios_weib50.csv') )

### Comparative simulations for 50th percentile

#-------------------- Prep


#### Constants

set.seed(79285)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
M = 12
n = 40
lostart = 3
histart = 10

# Generate the scenarios on the dose grid
weib50F = weib50parm[1:nsim, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib50Fu = weib50parm[1:nsim, ][ ,as.vector(pweib3((1:M)-3, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib50Fl = weib50parm[1:nsim, ][ ,as.vector(pweib3((1:M)+3, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh50w = matrix(runif(n*nsim), nrow=n)
					

	
#-------------------- simulation and estimation: Weibull	

########### k-row standard	

w50lomid = dfsim(n, starting = lostart, Fvals=weib50F, ensemble=nsim,
						thresholds = thresh50w)					
estw50lomid = estbatch(w50lomid, truth=weib50parm$t50[1:nsim], target=.5)
	
# Good pausing point to see the output of one framework 
# stop('ha')	

w50himid = dfsim(n, starting = histart, Fvals=weib50F, ensemble=nsim,
						thresholds = thresh50w)					
estw50himid = estbatch(w50himid, truth=weib50parm$t50[1:nsim], target=.5)

w50midmid = dfsim(n, starting = M/2, Fvals=weib50F, ensemble=nsim,
						thresholds = thresh50w)					
estw50midmid = estbatch(w50midmid, truth=weib50parm$t50[1:nsim], target=.5)

w50minhi = dfsim(n, starting = 1, Fvals=weib50Fu, ensemble=nsim,
						thresholds = thresh50w)					
estw50minhi = estbatch(w50minhi, truth=weib50parm$t50[1:nsim]+3, target=.5)
	
w50maxlo = dfsim(n, starting = M, Fvals=weib50Fl, ensemble=nsim,
						thresholds = thresh50w)					
estw50maxlo = estbatch(w50maxlo, truth=weib50parm$t50[1:nsim]-3, target=.5)
				

save.image(file.path(outdir, 'grandsim50w.RData'))					
cat('done.\n')					
									





		

		
					





