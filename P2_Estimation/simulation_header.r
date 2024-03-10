# rm(list=ls())
library(cir)
library(upndown)
library(data.table)
setDTthreads(2)

####----------------------- Constants and quick Utilities

outdir = '../../output'
# source('code_for_Assaf.R')

#### Simple performance metrics, with some twists

rmse = function(x,ref,na.rm = TRUE, winsor = FALSE) 
{
	tmp = (x-ref)^2
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

####----------------------- Functions to generate the dose-response scenarios

weibshift <- function(shp, scl, targx = 5.5, targy)
{
	qweibull(targy, shape = shp, scale = scl) - targx
}
pweib3 <- function(x, shp, scl, shift) pweibull(q=x+shift, shape=shp, scale=scl)
qweib3 <- function(p, shp, scl, shift) qweibull(p, shape=shp, scale=scl) - shift


####---------------- Batch calculation of estimators and their metrics  ###########

### Parallelized for Windows environment via 'foreach'

estbatch <- function(simdat, truth, target, bpt=target, rawout=TRUE, cores = 3, 
			n = NULL, nsim = NULL,  B = 250, randboot = TRUE, 
			cirb = FALSE, desfun=krow, desargs=list(k=1), 
			doseset = NULL, conf = 0.9, bigerr = 0.9, addLiao = FALSE)

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
if(is.null(doseset)) doseset = 1:M
if(length(doseset) != M) stop('Mistmatch in length of dose set.\n')
# interval tails
ctail = (1 - conf) / 2

# If all identical and truth is a scalar:
if(length(truth) == 1) truth = rep(truth, nsim)
truth = truth[1:nsim]
if(length(truth) != nsim) stop('Mistmatch in length of true values.\n')


ests <- foreach(a = 1:nsim, .combine = 'rbind', 
			.packages = c('cir','upndown','plyr') ) %dopar%  {
### First, generating the bootstrap sample for all CI estimation
#      (and we get dynamean() bootstrap CI's "for free")
	eout = data.frame(true = truth[a])
	boots = dfboot(simdat$doses[1:n, a], simdat$responses[1:n,a], B=B, doses = doseset,
				design = desfun, desArgs = desargs, showdots = FALSE,
				target = target, balancePt = bpt, full = TRUE, randstart = randboot)
#				return(boots)

# "Dressing up" the dose levels (which are 1:m in the progress loop above) with real values
    if(any(doseset != 1:M)) bootdoses = suppressMessages(plyr::mapvalues(boots$x, 1:M, doseset) ) else bootdoses = boots$x

### Averaging estimators
	eout$dm48 = dixonmood(simdat$doses[1:n, a], simdat$response[1:n,a])
	eout$all1 = reversmean(simdat$doses[1:(n+1),a],simdat$response[1:n,a],rstart=1, conf = NULL)
	eout$all3 = reversmean(simdat$doses[1:(n+1),a],simdat$response[1:n,a],rstart=3, conf = NULL)
# Wetherill's estimator
	eout$rev1 = reversmean(simdat$doses[1:n,a],simdat$response[1:n,a],rstart=1, all=FALSE, conf = NULL)
	eout$dyna = dynamean(simdat$doses[1:(n+1),a], maxExclude = 1/2, conf = NULL)

### isotonics
	tmp1 = try(udest(simdat$doses[1:n, a], simdat$response[1:n,a], target=target, 
	balancePt = bpt, conf = conf) )
	tmp2 = try(udest(simdat$doses[1:n, a], simdat$response[1:n,a], target=target, 
	balancePt = bpt, conf = conf, estfun = oldPAVA) )
	eout$cir = ifelse(is.finite(tmp1$point), tmp1$point, NA)
	eout$ir = ifelse(is.finite(tmp2$point),tmp2$point, NA)

#### CIs

# Isotonic is simplest:
	eout$cirl = ifelse(is.finite(tmp1$point), unlist(tmp1[3]), NA)
	eout$ciru = ifelse(is.finite(tmp1$point), unlist(tmp1[4]), NA)
	eout$irl = ifelse(is.finite(tmp2$point), unlist(tmp2[3]), NA)
	eout$iru = ifelse(is.finite(tmp2$point), unlist(tmp2[4]), NA)
	if(addLiao)
	{
		tmp3 = udest(simdat$doses[1:n, a], simdat$response[1:n,a], target=target, 
				balancePt = bpt, conf = conf, intfun = liaoCI)
		eout$liaol = unlist(tmp3[3])
		eout$laiuu = unlist(tmp4[4])
	}

# Dynamic mean also already available via the dfboot call:
	
	tmp = quantile(boots$ests, probs = c(ctail, 1-ctail), type = 6)
	eout$dynal = tmp[1]
	eout$dynau = tmp[2]

# Bootstrap for the other estimators

	all3boot = rep(NA, B)
	rev1boot = all3boot
	cirboot = all3boot
	for(b in 1:B) 
	{
			all3boot[b] = reversmean(x = bootdoses[,b], y = boots$y[,b], 
                       conf = NULL)			
			rev1boot[b] = reversmean(x = bootdoses[,b], y = boots$y[,b], 
                        all = FALSE, rstart = 1, conf = NULL)						
			if(cirb) {
				cirboot[b] = try(udest(x = bootdoses[1:n,b], y = boots$y[,b], 
                        conf = NULL, target = target, balancePt = bpt)	)
						if(!is.finite(cirboot[b])) cirboot[b] = NA
			}
	}	

	tmp = quantile(all3boot, probs = c(ctail, 1-ctail), type = 6, na.rm = TRUE)
	eout$all3l = tmp[1]
	eout$all3u = tmp[2]
	tmp = quantile(rev1boot, probs = c(ctail, 1-ctail), type = 6, na.rm = TRUE)
	eout$rev1l = tmp[1]
	eout$rev1u = tmp[2]
	if(cirb) {
	tmp = quantile(cirboot, probs = c(ctail, 1-ctail), type = 6, na.rm = TRUE)
	eout$cbootl = tmp[1]
	eout$cbootu = tmp[2]
	}
	
	eout
}
stopCluster(cl)
cat(base::date(), '\n')

setDT(ests)

##### returning everything raw - now it's the default:

if(rawout) return(ests)

#### Otherwise returning headline summaries

### Point estimate performance
tmp=list(metrics=ests[ ,apply(.SD, 2, trio, ref=true, p=bigerr), .SDcol=names(ests)[2:8] ],
	cirmissed=mean(is.na(ests$cir)), n = n, ensemble = nsim, 
	target = target, startpt = simdat$doses[1,1], targloc = mean(truth))

### CI coverage

tmp$coverage = ests[ , list(all3 = mean(all3l<=true & all3u>=true, na.rm=TRUE),
	rev1 = mean(rev1l<=true & rev1u>=true, na.rm=TRUE),
	dyna = mean(dynal<=true & dynau>=true, na.rm=TRUE),
	cir = mean(cirl<=true & ciru>=true, na.rm=TRUE),
	ir = mean(irl<=true & iru>=true, na.rm=TRUE),
	cirboot = ifelse(cirb, NA, mean(cbootl<=true & cbootu>=true, na.rm=TRUE) )
	) ]

### CI width
	ests[ , cirfin := (is.finite(ciru) & is.finite(cirl) ) ]
	ests[ , irfin :=  (is.finite(iru) & is.finite(irl) ) ]
	
tmp$widths = ests[ , list(all3 = mean(all3u - all3l, na.rm=TRUE),
	rev1 = mean(rev1u - rev1l, na.rm=TRUE),
	dyna = mean(dynau - dynal, na.rm=TRUE),
	cir = mean(ciru[cirfin] - cirl[cirfin]),
	ir = mean(iru[irfin] - irl[irfin]),
	cfinite = mean(cirfin) , ifinite = mean(irfin),
#	liao = ifelse(addLiao, mean(liaou - liaol, na.rm=TRUE), NA),
	cirboot = ifelse(cirb, NA, mean(cbootu - cbootl, na.rm=TRUE) ) 
	) ] 

tmp
}


