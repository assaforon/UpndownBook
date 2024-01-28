rm(list=ls())
library(cir)
library(upndown)
library(data.table)

####----------------------- Constants and quick Utilities

outdir = '../../output'

#### Simple performance metrics
rmse = function(x,ref,na.rm=TRUE) sqrt(mean((x-ref)^2,na.rm=na.rm))
bias = function(x,ref,na.rm=TRUE) mean(x-ref,na.rm=na.rm)
# Quantile of absolute error
qae = function(x,ref,na.rm=TRUE, p = 0.95) quantile(abs(x-ref), probs = p, na.rm=na.rm)



####----------------------- Functions to generate the dose-response scenarios

weibshift <- function(shp, scl, targx = 5.5, targy)
{
	qweibull(targy, shape = shp, scale = scl) - targx
}
pweib3 <- function(x, shp, scl, shift) pweibull(q=x+shift, shape=shp, scale=scl)
qweib3 <- function(p, shp, scl, shift) qweibull(p, shape=shp, scale=scl) - shift


####---------------- Batch calculation of estimators and their metrics  ###########

estbatch <- function(simdat, truth, target, bpt=target, rawout=FALSE, dots = TRUE, 
               B = 10, desfun, desargs, doseset, conf = 0.9)

{
require(cir)
require(upndown)
require(plyr)
require(data.table)

sizes=dim(simdat$response)
n=sizes[1]
nsim=sizes[2]
M = dim(simdat$scenarios)[1]
if(length(doseset) != M) stop('Mistmatch in length of dose set.\n')
# interval tails
ctail = (1 - conf) / 2

ests = data.table(true=truth)
cis = copy(ests)

for (a in 1:nsim)
{
### First, generating the bootstrap sample for all CI estimation
#      (and we get dynamean() bootstrap "for free")

	boots = dfboot(simdat$dose[1:n, a], simdat$response[ ,a], B=B, doses = doseset,
				design = desfun, desArgs = desargs, target = target, balancePt = bpt)
# "Dressing up" the dose levels (which are 1:m in the progress loop above) with real values
    bootdoses = suppressMessages(plyr::mapvalues(boots$doses, 1:M, doseset) )

### Averaging estimators
	ests$dm48[a] = dixonmood(simdat$dose[1:n, a], simdat$response[ ,a])
	ests$all1[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=1, conf = NULL)
	ests$all3[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=3, conf = NULL)
# Wetherill's estimator
	ests$rev1[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=1, , conf = NULL)
	ests$dyna[a] = dynamean(simdat$dose[,a], maxExclude = 1/2, conf = NULL)
	
### isotonics
	tmp1 = udest(simdat$dose[1:n, a], simdat$response[ ,a], target=target, 
	balancePt = bpt, conf = conf)
	tmp2 = udest(simdat$dose[1:n, a], simdat$response[ ,a], target=target, 
	balancePt = bpt, conf = conf, estfun = oldPAVA)
	ests$cir[a] = tmp1$point
	ests$ir[a] = tmp2$point

#### CI

# Isotonic is simplest:
	cis$cirl[a] = tmp1[3]
	cis$ciru[a] = tmp1[4]
	cis$irl[a] = tmp2[3]
	cis$iru[a] = tmp2[4]

	tmp = quantile(boots$ests, probs = c(tail, 1-tail), type = 6)
	cis$dynal[a] = tmp[1]
	cis$dynau[a] = tmp[2]

# Bootstrap for the other avging
	all3boot = rep(NA, B)
	rev1boot = all3boot
	cirboot = all3boot
	for(b in 1:B) {
			all3boot[b] = reversmean(x = bootdoses[,b], y = bootdat$responses[,b], 
                       conf = NULL)			
			rev1boot[b] = reversmean(x = bootdoses[,b], y = bootdat$responses[,b], 
                        all = FALSE, rstart = 1, conf = NULL)						
			cirboot[b] = udest(x = bootdoses[1:n,b], y = bootdat$responses[,b], 
                        conf = NULL, target = target, balancePt = bpt)	
	}	

	tmp = quantile(all3boot, probs = c(tail, 1-tail), type = 6, na.rm = TRUE)
	cis$all3l[a] = tmp[1]
	cis$all3u[a] = tmp[2]
	tmp = quantile(rev1boot, probs = c(tail, 1-tail), type = 6, na.rm = TRUE)
	cis$rev1l[a] = tmp[1]
	cis$rev1u[a] = tmp[2]
	tmp = quantile(cirboot, probs = c(tail, 1-tail), type = 6, na.rm = TRUE)
	cis$cbootl[a] = tmp[1]
	cis$cbootu[a] = tmp[2]

	if(dots & a%%10==0) cat('.')
	if(dots & a%%100==0) cat('\n')
}


# returning everything:
if(rawout) return(list(point=ests, ci=cis))
# returning headline summaries:
tmp=list(metrics=ests[,apply(.SD,2,duo,ref=true,na.rm=TRUE),.SDcol=names(ests)[-1]],
	irmissed=mean(is.na(ests$ir)))
if(ci) ## CI performance
{
	tmp$coverage=cis[,list(all3n=mean(all3lo0<=true & all3hi0>=true,na.rm=TRUE),
	all3=mean(all3lo<=true & all3hi>=true,na.rm=TRUE),
	all3Wid=mean(all3hi-all3lo,na.rm=TRUE),
	cir=mean(cirl<=true & ciru>=true,na.rm=TRUE),
	ir=mean(irl<=true & iru>=true,na.rm=TRUE),
	cirWid=mean(ciru-cirl,na.rm=TRUE), cirWid_old=mean(ciru_old-cirl_old,na.rm=TRUE), 
	irWid=mean(iru-irl,na.rm=TRUE))]
}
tmp
}


