# rm(list=ls())
library(cir)
library(upndown)
library(data.table)

####----------------------- Constants and quick Utilities

outdir = '../../output'
# source('code_for_Assaf.R')

### New! Wetherill
# source('wetherill65.r')

#### Simple performance metrics, with some twists

rmse = function(x,ref,na.rm = TRUE, exclude = NULL, winsor = FALSE) 
{
	tmp = (x-ref)^2
	if(!is.null(exclude)) tmp = tmp[-exclude]
	if(winsor) tmp[is.na(tmp)] = max(tmp, na.rm = TRUE)
	sqrt(mean(tmp, na.rm=na.rm))
}

bias = function(x,ref,na.rm=TRUE) mean(x-ref,na.rm=na.rm)
# Trimmed mean absolute error 
mae = function(x,ref,na.rm=TRUE, p = 0.9, winsor = TRUE) 
{
	n = length(x)
	tmp = sort( abs(x-ref) , na.last = TRUE)
	if(winsor) tmp[is.na(tmp)] = max(tmp, na.rm = TRUE)
 mean(tmp[1:round(n*p)], na.rm = na.rm) 
}
# Combos galore!
duo = function(x,ref,na.rm=TRUE) c(rmse=rmse(x,ref), bias=bias(x,ref) )
trio = function(x,ref,na.rm=TRUE, p=0.9, ...) {
	tmp=c(rmse=rmse(x,ref, ...), bias=bias(x,ref), mae = mae(x,ref, p=p, ...) )
	names(tmp) = c('RMSE', 'Bias', paste('MAE', round(p*100), sep='') )
	tmp
}

###--------------------------- Phase-I specific utilities

### MTD selection
#     Function returns an index of the "MTD" estimated dose
# Fests: vector of F-hat estimates via the design's estimator of choice
# To get the largest dose <= target, set exclude=target.
whichmtd<-function(Fests , targ, cutoff=1,...) 
{
#print(Fests)
if(min(Fests,na.rm=TRUE)>=cutoff) {
#	cat("worked?")
	return(0) 
 } else return(which.min(abs(Fests[Fests<cutoff]-targ)))
}


### Pre-filtering utility for scenarios
inwindow=function(Fvals,lo,hi)
{
	if (lo>=hi) stop('low bound not lower than high bound!')
	return(sum(Fvals >= lo & Fvals <= hi) )
}

## How many treated on MTD?

matnstar=function(allocs,Fmat, target) mapply(function(x,y) {sum(x==y & !is.na(x))}, x = split(allocs,col(allocs)),
	 y = apply(Fmat, 2, function(x,targ) which.min(abs(x-targ)), targ=target) )

## How many treated inside interval?
matninterval=function(allocs,Fmat,lo,hi) mapply(function(x,y,lo,hi) 
{refs=which(y>=lo & y<=hi);sum(x %in% refs & !is.na(x))},split(allocs,col(allocs)),split(Fmat,col(Fmat)),MoreArgs=list(lo=lo,hi=hi))

estCRM<-function(doses, responses, target=0.3, cutoff=1, getF = FALSE, pointest = TRUE, skel)
{
require(dfcrm)
y=crm(prior=skel, target=target, level=doses, tox=responses)$ptox
if(getF) return(y)
if(pointest) pest = approx(y, 1:length(skel), xout = target)$y 

include = 1:length(y)
if(is.finite(cutoff)) include = which(y<=cutoff)

mtd = which.min(abs(y-target))
if(length(mtd) == 0) mtd = NA

if(pointest) return(c(point=pest, mtd=mtd))
return(mtd)
}

# stop('ya')
####----------------------- Functions to generate the dose-response scenarios

weibshift <- function(shp, scl, targx = 5.5, targy)
{
	qweibull(targy, shape = shp, scale = scl) - targx
}
pweib3 <- function(x, shp, scl, shift) pweibull(q=x+shift, shape=shp, scale=scl)
qweib3 <- function(p, shp, scl, shift) qweibull(p, shape=shp, scale=scl) - shift

####---------------- Batch calculation of estimators and their metrics  ###########
### Compatible with UDDs and interval designs
### Parallelized for Windows environment via 'foreach'

restbatch <- function(simdat, truth, target, bpt=target, halfwidth, cores = 13, n = NULL, nsim = NULL, conf = 0.9)
{
cat(base::date(), '\n')
require(cir)
require(upndown)
require(plyr)
require(data.table)
require(doParallel)
cl <- makeCluster(cores, type = "SOCK")
registerDoParallel(cl)

sizes=dim(simdat$response)
if(is.null(n)) n=sizes[1]
if(is.null(nsim)) nsim=sizes[2]
M = dim(simdat$scenarios)[1]

# If all identical and truth is a scalar:
if(length(truth) == 1) truth = rep(truth, nsim)
truth = truth[1:nsim]
if(length(truth) != nsim) stop('Mistmatch in length of true values.\n')

#--------------------- Parallel loop
ests0 <- foreach(a = 1:nsim,  
			.packages = c('cir','upndown','plyr') ) %dopar%  {
	eout = data.frame(true = truth[a])
	
### Number within a tolerance interval
	goodF = which(simdat$scenarios[ ,a] >= target-halfwidth & simdat$scenarios[ ,a] <= target+halfwidth)
	eout$ninterval = sum(!is.na(simdat$doses[1:n, a]) & simdat$doses[1:n, a] %in% goodF )
	
	if(var(simdat$doses[1:n, a]) == 0) return(eout) # a dud run, nothing else can be calculated

### isotonics
	tmp1 = try(udest(simdat$doses[1:n, a], simdat$responses[1:n,a], target=target, 
	balancePt = bpt, conf = conf) )

	eout$pointest = ifelse( 'data.frame' %in% class(tmp1), tmp1$point, NA)

### Using isotonics for ordinal MTD estimate
	tmp2 = try(quickIsotone(x=simdat$doses[1:n, a], y=simdat$responses[1:n,a], target=bpt, 
	conf = conf, adaptiveShrink = TRUE) )

	eout$mtdest = ifelse( 'data.frame' %in% class(tmp2), tmp2$x[which.min(abs(tmp2$y-target))], NA)
### Est within a tolerance interval
	eout$goodest = (eout$mtdest %in% goodF)
		
	eout	
}
stopCluster(cl)
ests = rbindlist(ests0, fill = TRUE) # Also converts to data.table

ests[ , mtd := apply(simdat$scenarios, 2, function(x,targ) which.min(abs(x-targ)), targ=target) ]
tmp = simdat$doses[1:n, ]
ests[ , nstar := mapply(function(x,y) {sum(x==y & !is.na(x))}, x = split(tmp,col(tmp)), y = mtd) ]
cat(base::date(), '\n')
return(ests)
}


####---------------- Batch calculation of CRM estimators and their metrics  ###########
### Parallelized for Windows environment via 'foreach'

crmbatch <- function(simdat, truth, target, skel, halfwidth, cores = 13, n = NULL, nsim = NULL,  conf = 0.9)
{
cat(base::date(), '\n')
require(dfcrm)
require(plyr)
require(data.table)
require(doParallel)
cl <- makeCluster(cores, type = "SOCK")
registerDoParallel(cl)

sizes=dim(simdat$response)
if(is.null(n)) n=sizes[1]
if(is.null(nsim)) nsim=sizes[2]
M = dim(simdat$scenarios)[1]

# If all identical and truth is a scalar:
if(length(truth) == 1) truth = rep(truth, nsim)
truth = truth[1:nsim]
if(length(truth) != nsim) stop('Mistmatch in length of true values.\n')

#--------------------- Parallel loop
ests0 <- foreach(a = 1:nsim,  
			.packages = c('dfcrm','plyr') ) %dopar%  {
	eout = data.frame(true = truth[a])
	
### Number within a tolerance interval
	goodF = which(simdat$scenarios[ ,a] >= target-halfwidth & simdat$scenarios[ ,a] <= target+halfwidth)
	eout$ninterval = sum(!is.na(simdat$doses[1:n, a]) & simdat$doses[1:n, a] %in% goodF )

	if(var(simdat$doses[1:n, a]) == 0) return(eout) # a dud run, will trigger errors if attempted

	tmp = crm(prior=skel, target=target, level=simdat$doses[1:n, a], tox=simdat$responses[ , a])
	eout$pointest = approx(tmp$ptox, 1:M, xout = target)$y 
	eout$mtdest = which.min(abs(tmp$ptox - target) )
	eout$thetahat = tmp$estimate
### Est within a tolerance interval
	eout$goodest = (eout$mtdest %in% goodF)
		
	eout	
}
stopCluster(cl)
ests = rbindlist(ests0, fill = TRUE) # Also converts to data.table

ests[ , mtd := apply(simdat$scenarios, 2, function(x,targ) which.min(abs(x-targ)), targ=target) ]
tmp = simdat$doses[1:n, ]
ests[ , nstar := mapply(function(x,y) {sum(x==y & !is.na(x))}, x = split(tmp,col(tmp)), y = mtd) ]
cat(base::date(), '\n')
return(ests)
}



