cat(base::date(), '\n')
rm(list=ls())
library(magrittr)
library(upndown)
library(dfcrm)
library(BOIN)
library(data.table)

#---------------------- Preamble

source('OtherDesigns.r')
source('simother_header.r')
outdir = 'C:/GitHub/output'

weib90parm = fread(file.path(outdir, 'scenarios_weib90.csv') )

### Comparative simulations for 90th percentile

#-------------------- Prep

#### Constants

set.seed(79286)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
overhead = 2 # additional candidate scenarios, for some more filtering
M = 12
n = 50
targ0 = 0.9
goodint = 0.075 # half-width of a "good" interval around target

### Pre-filtering utility for scenarios
mindens = 0
maxdens = 4
tootight = 0.1
toolax = 0.3

ktarg90 = k2targ(6)
k90list = list(k=6, lowTarget=FALSE, fastStart=TRUE)
lostart = 3
histart = 10
#b90list = list(coin=3/7, lowTarget=TRUE, fastStart=FALSE)
### Info gathering for the other designs	
# Various CRM prior skeletons		
skel10_025 = getprior(halfwidth = 0.025, target = targ0, nu = 10, nlevel = M)
skel03_025 = getprior(halfwidth = 0.025, target = targ0, nu = 3, nlevel = M)
skel07_05 = getprior(halfwidth = 0.05, target = targ0, nu = 7, nlevel = M)
skel07_025 = getprior(halfwidth = 0.025, target = targ0, nu = 7, nlevel = M)

ccd90list = list(hwidth=0.09, targ=targ0)  # As per article!
ccd90list2 = list(hwidth=0.05, targ=targ0)


# Generate the scenarios on the dose grid

weib90F0 = weib90parm[1:(overhead*nsim), ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib90Fu0 = weib90parm[1:(overhead*nsim), ][ ,as.vector(pweib3((1:M)-3, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)

tightness = weib90F0[8, ] - weib90F0[5, ]
straddle = weib90F0[7, ] - weib90F0[6, ]

Fdensities = apply(weib90F0, 2, inwindow, lo=targ0-goodint, hi=targ0+goodint)
Fudensities = apply(weib90Fu0, 2, inwindow, lo=targ0-goodint, hi=targ0+goodint)

scen0 = which(Fdensities %in% mindens:maxdens & tightness > tootight & straddle < toolax)
scenu = which(Fudensities %in% mindens:maxdens & tightness > tootight & straddle < toolax)

weib90F = weib90F0[ ,scen0][ , 1:nsim]
weib90Fu = weib90Fu0[ , scenu][ , 1:nsim]
# For the record (and same for mid and hi targets)
densities = apply(weib90F, 2, inwindow, lo=targ0-goodint, hi=targ0+goodint)

truew90 = weib90parm[1:(overhead*nsim), ][scen0, ][1:nsim, t90bar]
truew90u = weib90parm[1:(overhead*nsim), ][scenu, ][1:nsim, t90bar] + 3



# stop('prep')
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh90w = matrix(runif(n*nsim), nrow=n)
			
#-------------------- simulation and estimation: Weibull	

kw90himid = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
restkw90himid = restbatch(kw90himid, truth=truew90, target=targ0, bpt=ktarg90, halfwidth = goodint)

ccdw90himid = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list, thresholds = thresh90w)					
restccdw90himid = restbatch(ccdw90himid, truth=truew90, target=targ0, halfwidth = goodint)
ccdw90himid2 = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list2, thresholds = thresh90w)					
restccdw90himid2 = restbatch(ccdw90himid2, truth=truew90, target=targ0, halfwidth = goodint)
				
crmw90himid10_025 = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel10_025, targ = targ0), thresholds = thresh90w)					
restcrmw90himid10_025 = crmbatch(crmw90himid10_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel10_025)
crmw90himid03_025 = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel03_025, targ = targ0), thresholds = thresh90w)					
restcrmw90himid03_025 = crmbatch(crmw90himid03_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel03_025)
crmw90himid07_05 = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_05, targ = targ0), thresholds = thresh90w)					
restcrmw90himid07_05 = crmbatch(crmw90himid07_05, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_05)
crmw90himid07_025 = dfsim(n, starting = histart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_025, targ = targ0), thresholds = thresh90w)					
restcrmw90himid07_025 = crmbatch(crmw90himid07_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_025)

cat('himid\n')	
# Good pausing point to see the output of one framework 
 
# stop('ha')	

kw90midmid = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
restkw90midmid = restbatch(kw90midmid, truth=truew90, target=targ0, bpt=ktarg90, halfwidth = goodint)

ccdw90midmid = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list, thresholds = thresh90w)					
restccdw90midmid = restbatch(ccdw90midmid, truth=truew90, target=targ0, halfwidth = goodint)
ccdw90midmid2 = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list2, thresholds = thresh90w)					
restccdw90midmid2 = restbatch(ccdw90midmid2, truth=truew90, target=targ0, halfwidth = goodint)

crmw90midmid10_025 = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel10_025, targ = targ0), thresholds = thresh90w)					
restcrmw90midmid10_025 = crmbatch(crmw90midmid10_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel10_025)
crmw90midmid03_025 = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel03_025, targ = targ0), thresholds = thresh90w)					
restcrmw90midmid03_025 = crmbatch(crmw90midmid03_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel03_025)
crmw90midmid07_05 = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_05, targ = targ0), thresholds = thresh90w)					
restcrmw90midmid07_05 = crmbatch(crmw90midmid07_05, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_05)
crmw90midmid07_025 = dfsim(n, starting = M/2, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_025, targ = targ0), thresholds = thresh90w)					
restcrmw90midmid07_025 = crmbatch(crmw90midmid07_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_025)

cat('midmid\n')
	
kw90minhi = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
restkw90minhi = restbatch(kw90minhi, truth=truew90u, target=targ0, bpt=ktarg90, halfwidth = goodint)

ccdw90minhi = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design = ccd,
						desArgs=ccd90list, thresholds = thresh90w)					
restccdw90minhi = restbatch(ccdw90minhi, truth=truew90u, target=targ0, halfwidth = goodint)
ccdw90minhi2 = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design = ccd,
						desArgs=ccd90list2, thresholds = thresh90w)					
restccdw90minhi2 = restbatch(ccdw90minhi2, truth=truew90u, target=targ0, halfwidth = goodint)
				
crmw90minhi10_025 = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel10_025, targ = targ0), thresholds = thresh90w)					
restcrmw90minhi10_025 = crmbatch(crmw90minhi10_025, truth=truew90u, target=targ0, halfwidth = goodint, skel = skel10_025)
crmw90minhi03_025 = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel03_025, targ = targ0), thresholds = thresh90w)					
restcrmw90minhi03_025 = crmbatch(crmw90minhi03_025, truth=truew90u, target=targ0, halfwidth = goodint, skel = skel03_025)
crmw90minhi07_05 = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_05, targ = targ0), thresholds = thresh90w)					
restcrmw90minhi07_05 = crmbatch(crmw90minhi07_05, truth=truew90u, target=targ0, halfwidth = goodint, skel = skel07_05)
crmw90minhi07_025 = dfsim(n, starting = 1, Fvals=weib90Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_025, targ = targ0), thresholds = thresh90w)					
restcrmw90minhi07_025 = crmbatch(crmw90minhi07_025, truth=truew90u, target=targ0, halfwidth = goodint, skel = skel07_025)

cat('minhi\n')	

kw90lomid = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim,
						desArgs=k90list, thresholds = thresh90w)					
restkw90lomid = restbatch(kw90lomid, truth=truew90, target=targ0, bpt=ktarg90, halfwidth = goodint)

ccdw90lomid = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list, thresholds = thresh90w)					
restccdw90lomid = restbatch(ccdw90lomid, truth=truew90, target=targ0, halfwidth = goodint)
ccdw90lomid2 = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design = ccd,
						desArgs=ccd90list2, thresholds = thresh90w)					
restccdw90lomid2 = restbatch(ccdw90lomid2, truth=truew90, target=targ0, halfwidth = goodint)

crmw90lomid10_025 = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel10_025, targ = targ0), thresholds = thresh90w)					
restcrmw90lomid10_025 = crmbatch(crmw90lomid10_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel10_025)
crmw90lomid03_025 = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel03_025, targ = targ0), thresholds = thresh90w)					
restcrmw90lomid03_025 = crmbatch(crmw90lomid03_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel03_025)
crmw90lomid07_05 = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_05, targ = targ0), thresholds = thresh90w)					
restcrmw90lomid07_05 = crmbatch(crmw90lomid07_05, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_05)
crmw90lomid07_025 = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel07_025, targ = targ0), thresholds = thresh90w)					
restcrmw90lomid07_025 = crmbatch(crmw90lomid07_025, truth=truew90, target=targ0, halfwidth = goodint, skel = skel07_025)

save.image(file.path(outdir, 'othersim90w.RData'))					
cat('Done.\n')					
									





		

		
					





