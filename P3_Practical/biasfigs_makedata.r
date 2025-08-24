rm(list=ls())
library(cir)
library(upndown)
library(dfcrm)
library(BOIN)
library(magrittr)
library(data.table)

# For use on your own machine, change 'outdir' to where you want the figures saved
outdir = ('../../output')

#------------ Prep

load(file.path(outdir, 'grandsim50l.RData'))	

n = 40
n2 = 100
nsim = 1000

logi50F = logi50parm[1:nsim, ][ ,as.vector(plogis(1:M, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
thresh50l = matrix(runif(n*nsim), nrow=n)
thresh50l_long = matrix(runif(n*nsim), nrow=n2)

## prepping for other designs
source('OtherDesigns.r')
skel6_05 = getprior(halfwidth = 0.05, target = 0.5, nu = M/2, nlevel = M)

boin10list = list(lookup = get.boundary(target=0.5, ncohort=n, cohortsize=1, cutoff.eli=1)$boundary_tab)
boin10list_long = list(lookup = get.boundary(target=0.5, ncohort=n2, cohortsize=1, cutoff.eli=1)$boundary_tab)

### Utility for one-stop-shop bias calculation

fcalc <- function(simdat, directbias = TRUE)
{
	m = dim(simdat$scenarios)[1]
	n = dim(simdat$responses)[1]
#	dcuts = cut(simdat$doses[1:n, ], (0:m)-0.5, labels = 1:m)
	tmp = mapply(function(y,x) sapply(split(y,x), mean), x=split(simdat$doses[1:n,], col(simdat$doses[1:n,])), 
							y=split(simdat$responses, col(simdat$responses)) )
							
	dout = sapply(tmp, function(x,mm) {dout = rep(NA, mm); dout[as.integer(names(x))] <- x; dout}, mm=m)
	if(directbias) return(rowMeans(dout - simdat$scenarios, na.rm=TRUE) )
	dout
}

#------------- Redo UD with B=1000, and add non-UD simulations

l50midmid = dfsim(n, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l)					

crml50midmid = dfsim(n, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l, design=wrapCRM,
						desArgs=list(skel = skel6_05, targ = 0.5))					

boinl50midmid = dfsim(n, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l, design=boin,
						desArgs=boin10list)	

l50midmid_long = dfsim(n2, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l_long)					

crml50midmid_long = dfsim(n2, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l_long, design=wrapCRM,
						desArgs=list(skel = skel6_05, targ = 0.5))					

boinl50midmid_long = dfsim(n2, starting = M/2, Fvals=logi50F, ensemble=nsim,
						thresholds = thresh50l_long, design=boin,
						desArgs=boin10list_long)	



#crml50bias = data.table(Design = 'CRM', Dose = 1:M, Bias = fcalc(crml50midmid) )			
#boinl50bias = data.table(Design = 'Interval (BOIN)', Dose = 1:M, Bias = fcalc(boinl50midmid) )				
#udl50bias = data.table(Design = 'Classical UDD', Dose = 1:M, Bias = fcalc(l50midmid) )

#biasdat = rbind(udl50bias, crml50bias, boinl50bias)

save(l50midmid, crml50midmid, boinl50midmid, 
		l50midmid_long, crml50midmid_long, boinl50midmid_long, 
		file = file.path(outdir, 'biasfigs_makedata.r') ) 		

		
