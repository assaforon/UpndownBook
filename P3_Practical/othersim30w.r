cat(base::date(), '\n')
rm(list=ls())
library(magrittr)
library(upndown)
library(dfcrm)
library(data.table)

#---------------------- Preamble

source('C:/GitHub/UpndownBook/P2_Estimation/simulation_header.r')

# source('simother_header.r')
source('OtherDesigns.r')
outdir = 'C:/GitHub/output'

weib30parm = fread(file.path(outdir, 'scenarios_weib30.csv') )

### Comparative simulations for 30th percentile

#-------------------- Prep


#### Constants

set.seed(379286)
nsim = 10 # This overrides 'nsim' in the scenario-generation scripts
M = 8
n = 30
ktarg30 = k2targ(2, lowTarget=TRUE)
k30list = list(k=2, lowTarget=TRUE, fastStart=FALSE)
#b30list = list(coin=3/7, lowTarget=TRUE, fastStart=FALSE)
### Info gathering for the other designs	
# Various CRM prior skeletons		
skel4_10 = getprior(halfwidth = 0.1, target = 0.3, nu = 4, nlevel = M)
skel4_05 = getprior(halfwidth = 0.05, target = 0.3, nu = 4, nlevel = M)
skel2_05 = getprior(halfwidth = 0.05, target = 0.3, nu = 2, nlevel = M)
skel6_05 = getprior(halfwidth = 0.05, target = 0.3, nu = 6, nlevel = M)

boin10list = list(lookup = get.boundary(target=0.3,ncohort=30,cohortsize=1)$boundary_tab)
ccd10list = list(hwidth=0.1, targ=0.3)

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
# estkw30minmid = estbatch(kw30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
#					bpt=ktarg30, desfun=krow, ccurvy = TRUE, desargs=k30list)

crmw30minmid4_10 = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_10, targ = 0.3, pace = 2), thresholds = thresh30w)					
crmw30minmid4_05 = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_05, targ = 0.3, pace = 2), thresholds = thresh30w)					
crmw30minmid4_05fast = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_05, targ = 0.3), thresholds = thresh30w)					
crmw30minmid2_05 = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel2_05, targ = 0.3, pace = 2), thresholds = thresh30w)					
crmw30minmid6_05 = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel6_05, targ = 0.3, pace = 2), thresholds = thresh30w)					
boinw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design = boin,
						desArgs=boin10list, thresholds = thresh30w)					
ccdw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim, design = ccd,
						desArgs=ccd10list, thresholds = thresh30w)					
	
# Good pausing point to see the output of one framework 
 
stop('ha')	

kw30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30midmid = estbatch(kw30midmid, truth=weib30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=krow, ccurvy = TRUE, desargs=k30list)

kw30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30minhi = estbatch(kw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=krow, ccurvy = TRUE, desargs=k30list)
# estkw30minhi_2 = estbatch(kw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
	#				bpt=ktarg30, desfun=krow, ccurvy = TRUE, desargs=k30list, randboot = FALSE)

kw30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim,
						desArgs=k30list, thresholds = thresh30w)					
estkw30minlo = estbatch(kw30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=krow, ccurvy = TRUE, desargs=k30list)
	
cat('k standard\n')	
					
########### BCD standard	

bw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)					
estbw30minmid = estbatch(bw30minmid, truth=weib30parm$t30[1:nsim], target=.3, 
					 desfun=bcd, ccurvy = TRUE, desargs=b30list)
		
bw30midmid = dfsim(n, starting = M/2, Fvals=weib30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30midmid = estbatch(bw30midmid, truth=weib30parm$t30[1:nsim], target=.3,  
				desfun=bcd, ccurvy = TRUE, desargs=b30list)

bw30minhi = dfsim(n, starting = 1, Fvals=weib30Fu, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30minhi = estbatch(bw30minhi, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					 desfun=bcd, ccurvy = TRUE, desargs=b30list)

bw30minlo = dfsim(n, starting = 1, Fvals=weib30Fl, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)			
estbw30minlo = estbatch(bw30minlo, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					 desfun=bcd, ccurvy = TRUE, desargs=b30list)

cat('b standard\n')						

save.image(file.path(outdir, 'grandsim30w.RData'))					
cat('b variants and done.\n')					
									





		

		
					





