cat(base::date(), '\n')
rm(list=ls())
library(magrittr)
library(upndown)
library(dfcrm)
library(BOIN)
library(data.table)

#---------------------- Preamble

#source('C:/GitHub/UpndownBook/P2_Estimation/simulation_header.r')

source('OtherDesigns.r')
source('simother_header.r')
outdir = 'C:/GitHub/output'

weib30parm = fread(file.path(outdir, 'scenarios_weib30.csv') )

### Comparative simulations for 30th percentile

#-------------------- Prep

#### Constants

set.seed(379286)
nsim =  500 # This overrides 'nsim' in the scenario-generation scripts
overhead = 2 # additional candidate scenarios, for some more filtering
M = 8
n = 30
csize = 3 # cohort size
targ0 = 0.3

### Pre-filtering utility for scenarios

maxdens = 3
tootight = 0.15

### Info gathering for the other designs	
# Various CRM prior skeletons		
skel4_10 = getprior(halfwidth = 0.1, target = targ0, nu = 4, nlevel = M)
skel2_05 = getprior(halfwidth = 0.05, target = targ0, nu = 2, nlevel = M)
skel7_05 = getprior(halfwidth = 0.05, target = targ0, nu = 7, nlevel = M)

boin10list = list(lookup = get.boundary(target=targ0,ncohort=30,cohortsize=1)$boundary_tab)
ccd10list = list(hwidth=0.1, targ=targ0)

# Up-and-downs!
b30list = list(coin = 3/7, lowTarget = TRUE)
g30list = list(s = csize, ll = 0, ul = 2)
gtarg30 = g2targ(3,0,2)

lostart = 1

#-------------- Generate the scenarios on the dose grid

weib30F0 = weib30parm[1:(overhead*nsim), ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib30Fu0 = weib30parm[1:(overhead*nsim), ][ ,as.vector(pweib3((1:M)-2, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
weib30Fl0 = weib30parm[1:(overhead*nsim), ][ ,as.vector(pweib3((1:M)+2, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)

tightness = pmax(weib30F0[6, ] - weib30F0[4, ], weib30F0[5, ] - weib30F0[3, ])

Fdensities = apply(weib30F0, 2, inwindow, lo=targ0-0.1, hi=targ0+0.1)
Fudensities = apply(weib30Fu0, 2, inwindow, lo=targ0-0.1, hi=targ0+0.1)
Fldensities = apply(weib30Fl0, 2, inwindow, lo=targ0-0.1, hi=targ0+0.1)

weib30F = weib30F0[ , Fdensities %in% 1:maxdens & tightness > tootight][ , 1:nsim]
weib30Fu = weib30Fu0[ , Fudensities %in% 1:maxdens & tightness > tootight][ , 1:nsim]
weib30Fl = weib30Fl0[ , Fldensities %in% 1:maxdens & tightness > tootight][ , 1:nsim]

truew30 = weib30parm[1:(overhead*nsim), ][Fdensities %in% 1:maxdens & tightness > tootight, ][1:nsim, t30bar]
truew30u = weib30parm[1:(overhead*nsim), ][Fudensities %in% 1:maxdens & tightness > tootight, ][1:nsim, t30bar] + 2
truew30l = weib30parm[1:(overhead*nsim), ][Fldensities %in% 1:maxdens & tightness > tootight, ][1:nsim, t30bar] - 2

# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh30w = matrix(runif(n*nsim), nrow=n)
			
#-------------------- simulation and estimation: Weibull	

bw30minmid = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design = bcd,
						desArgs=b30list, thresholds = thresh30w)					
restbw30minmid = restbatch(bw30minmid, truth=truew30, target=targ0, bpt=targ0, halfwidth = 0.1)

gw30minmid = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design = groupUD,
						desArgs=g30list, thresholds = thresh30w)					
restgw30minmid = restbatch(gw30minmid, truth=truew30, target=targ0, bpt=gtarg30, halfwidth = 0.1)

boinw30minmid = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design = boin,
						desArgs=boin10list, thresholds = thresh30w)					
restboinw30minmid = restbatch(boinw30minmid, truth=truew30, target=targ0, halfwidth = 0.1)

ccdw30minmid = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design = ccd,
						desArgs=ccd10list, thresholds = thresh30w)					
restccdw30minmid = restbatch(ccdw30minmid, truth=truew30, target=targ0, halfwidth = 0.1)

crmw30minmid4_10 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_10, targ = targ0), thresholds = thresh30w)	
restcrmw30minmid4_10 = crmbatch(crmw30minmid4_10, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel4_10)
						
crmw30minmid2_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel2_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minmid2_05 = crmbatch(crmw30minmid2_05, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel2_05)
crmw30minmid7_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel7_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minmid7_05 = crmbatch(crmw30minmid7_05, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel7_05)

cat('minmid\n')	
# Good pausing point to see the output of one framework 
 
# stop('ha')	

bw30midmid = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design = bcd,
						desArgs=b30list, thresholds = thresh30w)					
restbw30midmid = restbatch(bw30midmid, truth=truew30, target=targ0, bpt=targ0, halfwidth = 0.1)

gw30midmid = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design = groupUD,
						desArgs=g30list, thresholds = thresh30w)					
restgw30midmid = restbatch(gw30midmid, truth=truew30, target=targ0, bpt=gtarg30, halfwidth = 0.1)

boinw30midmid = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design = boin,
						desArgs=boin10list, thresholds = thresh30w)					
restboinw30midmid = restbatch(boinw30midmid, truth=truew30, target=targ0, halfwidth = 0.1)

ccdw30midmid = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design = ccd,
						desArgs=ccd10list, thresholds = thresh30w)					
restccdw30midmid = restbatch(ccdw30midmid, truth=truew30, target=targ0, halfwidth = 0.1)

crmw30midmid4_10 = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_10, targ = targ0), thresholds = thresh30w)	
restcrmw30midmid4_10 = crmbatch(crmw30midmid4_10, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel4_10)
						
crmw30midmid2_05 = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel2_05, targ = targ0), thresholds = thresh30w)					
restcrmw30midmid2_05 = crmbatch(crmw30midmid2_05, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel2_05)
crmw30midmid7_05 = dfsim(n, starting = M/2, cohort = csize, Fvals=weib30F, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel7_05, targ = targ0), thresholds = thresh30w)					
restcrmw30midmid7_05 = crmbatch(crmw30midmid7_05, truth=truew30, target=targ0, halfwidth = 0.1, skel = skel7_05)

cat('midmid\n')
	
bw30minhi = dfsim(n, starting = 1, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design = bcd,
						desArgs=b30list, thresholds = thresh30w)					
restbw30minhi = restbatch(bw30minhi, truth=truew30u, target=targ0, bpt=targ0, halfwidth = 0.1)

gw30minhi = dfsim(n, starting = 1, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design = groupUD,
						desArgs=g30list, thresholds = thresh30w)					
restgw30minhi = restbatch(gw30minhi, truth=truew30u, target=targ0, bpt=gtarg30, halfwidth = 0.1)

boinw30minhi = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design = boin,
						desArgs=boin10list, thresholds = thresh30w)					
restboinw30minhi = restbatch(boinw30minhi, truth=truew30u, target=targ0, halfwidth = 0.1)

ccdw30minhi = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design = ccd,
						desArgs=ccd10list, thresholds = thresh30w)					
restccdw30minhi = restbatch(ccdw30minhi, truth=truew30u, target=targ0, halfwidth = 0.1)

crmw30minhi4_10 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_10, targ = targ0), thresholds = thresh30w)	
restcrmw30minhi4_10 = crmbatch(crmw30minhi4_10, truth=truew30u, target=targ0, halfwidth = 0.1, skel = skel4_10)
						
crmw30minhi2_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel2_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minhi2_05 = crmbatch(crmw30minhi2_05, truth=truew30u, target=targ0, halfwidth = 0.1, skel = skel2_05)
crmw30minhi7_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fu, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel7_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minhi7_05 = crmbatch(crmw30minhi7_05, truth=truew30u, target=targ0, halfwidth = 0.1, skel = skel7_05)

cat('minhi\n')	

bw30minlo = dfsim(n, starting = 1, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design = bcd,
						desArgs=b30list, thresholds = thresh30w)					
restbw30minlo = restbatch(bw30minlo, truth=truew30l, target=targ0, bpt=targ0, halfwidth = 0.1)

gw30minlo = dfsim(n, starting = 1, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design = groupUD,
						desArgs=g30list, thresholds = thresh30w)					
restgw30minlo = restbatch(gw30minlo, truth=truew30l, target=targ0, bpt=gtarg30, halfwidth = 0.1)

boinw30minlo = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design = boin,
						desArgs=boin10list, thresholds = thresh30w)					
restboinw30minlo = restbatch(boinw30minlo, truth=truew30l, target=targ0, halfwidth = 0.1)

ccdw30minlo = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design = ccd,
						desArgs=ccd10list, thresholds = thresh30w)					
restccdw30minlo = restbatch(ccdw30minlo, truth=truew30l, target=targ0, halfwidth = 0.1)

crmw30minlo4_10 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel4_10, targ = targ0), thresholds = thresh30w)	
restcrmw30minlo4_10 = crmbatch(crmw30minlo4_10, truth=truew30l, target=targ0, halfwidth = 0.1, skel = skel4_10)
						
crmw30minlo2_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel2_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minlo2_05 = crmbatch(crmw30minlo2_05, truth=truew30l, target=targ0, halfwidth = 0.1, skel = skel2_05)
crmw30minlo7_05 = dfsim(n, starting = lostart, cohort = csize, Fvals=weib30Fl, ensemble=nsim, design=wrapCRM,
						desArgs=list(skel = skel7_05, targ = targ0), thresholds = thresh30w)					
restcrmw30minlo7_05 = crmbatch(crmw30minlo7_05, truth=truew30l, target=targ0, halfwidth = 0.1, skel = skel7_05)

save.image(file.path(outdir, 'cohort3othersim30w.RData'))					
cat('Done.\n')					
									





		

		
					





