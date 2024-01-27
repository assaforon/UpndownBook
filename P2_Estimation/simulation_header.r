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

### Averaging estimators
	ests$dm48[a] = dixonmood(simdat$dose[1:n, a], simdat$response[ ,a])
	ests$all1[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=1, conf = NULL)
	ests$all3[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=3, conf = NULL)
# Wetherill's estimator
	ests$rev1[a] = reversmean(simdat$dose[,a],simdat$response[,a],rstart=1, , conf = NULL)
	ests$dyna[a] = dynamean(simdat$dose[,a], maxExclude = 1/2, conf = NULL)
	
### isotonics
	tmp1 = udest((simdat$dose[1:n, a], simdat$response[ ,a]) target=target, balancePt = bpt)
	tmp2 = udest((simdat$dose[1:n, a], simdat$response[ ,a]) target=target, balancePt = bpt, estfun = oldPAVA)
	ests$cir[a] = tmp1$point
	ests$ir[a] = tmp2$point

	
	ncand = max(table(tmpx))-1
		cis$all3neff[a] = ifelse(ncand>1, ncand, 2)
		
		cis$all3sd0[a]=sd(tmpx)
		cis$all3sd[a]=diff(quantile(tmpx,c(.1,.9),type=6))/2
		
		cis$all3se0[a]=cis$all3sd0[a]/sqrt(cis$all3neff[a])
		cis$all3lo0[a]=ests$all3[a]+qt(0.05,df=cis$all3neff[a]-1)*cis$all3se0[a]
		cis$all3hi0[a]=ests$all3[a]+qt(0.95,df=cis$all3neff[a]-1)*cis$all3se0[a]
		cis$all3se[a]=cis$all3sd[a]/sqrt(cis$all3neff[a])
		cis$all3lo[a]=ests$all3[a]+qt(0.05,df=cis$all3neff[a]-1)*cis$all3se[a]
		cis$all3hi[a]=ests$all3[a]+qt(0.95,df=cis$all3neff[a]-1)*cis$all3se[a]
	}

	if(ci) ## IR/CIR CIs
	{
		cis$cirl[a]=ifelse(!exists('tmp1') | !is.finite(tmp1$point),
			NA,tmp1$lower90conf)	
		cis$ciru[a]=ifelse(!exists('tmp1') | !is.finite(tmp1$point),
			NA,tmp1$upper90conf)	
		cis$cirl_old[a]=ifelse(!exists('tmp1') | !is.finite(tmp1$point),
			NA,tmp1_old$lower90conf)	
		cis$ciru_old[a]=ifelse(!exists('tmp1') | !is.finite(tmp1$point),
			NA,tmp1_old$upper90conf)	
		cis$irl[a]=ifelse(!exists('tmp2') | !is.finite(tmp2$point),
			NA,tmp2$lower90conf)	
		cis$iru[a]=ifelse(!exists('tmp2') | !is.finite(tmp2$point),
			NA,tmp2$upper90conf)	
	}
	if(dots & a%%10==0) cat('.')
	if(dots & a%%100==0) cat('\n')
}

if(irpush) # pushing NAs to get boundary values
{
	ests[is.na(cir),cir:=ifelse(ymin>target,1,ifelse(ymax<target,M,NA))]
	ests[is.na(ir),ir:=ifelse(ymin>target,1,ifelse(ymax<target,M,NA))]
}
ests[,c('ymin','ymax'):=NULL]


# returning everything:
if(rawout) return(list(point=ests,ci=cis))
# returning headline summaries:
tmp=list(metrics=ests[,apply(.SD,2,duo,ref=true,na.rm=TRUE),.SDcol=names(ests)[-1]],
	irmissed=mean(is.na(ests$ir)))
if(ci) ## CI performance
{
	tmp$coverage=cis[,list(all3n=mean(all3lo0<=true & all3hi0>=true,na.rm=TRUE),
	all3=mean(all3lo<=true & all3hi>=true,na.rm=TRUE),
	all3nWid=mean(all3hi0-all3lo0,na.rm=TRUE),all3Wid=mean(all3hi-all3lo,na.rm=TRUE),
	cir=mean(cirl<=true & ciru>=true,na.rm=TRUE),
	cir_old=mean(cirl_old<=true & ciru_old>=true,na.rm=TRUE),
	ir=mean(irl<=true & iru>=true,na.rm=TRUE),
	cirWid=mean(ciru-cirl,na.rm=TRUE), cirWid_old=mean(ciru_old-cirl_old,na.rm=TRUE), 
	irWid=mean(iru-irl,na.rm=TRUE))]
}
tmp
}


