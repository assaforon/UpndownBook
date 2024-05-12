####################################### Wrapper functions for the various designs  ###############################
## Compatible with upndown::dfsim()

### wrapping Cheung's dfcrm::crm() function to right argument order, with Goodman et al. restriction
wrapCRM <- function(doses, responses, skel, targ, pace=1, goodm=TRUE, goodmDown=TRUE, ...)
# Doses, responses: vectors with x and y values, respectively (y must be 0/1)
# skel: the F-curve skeleton to be raised to a power by the model. Best obtained ahead of time via dfcrm::getprior() 
# targ: target toxicity rate
# pace: the number of observations at same dose required before escalating. Default 1.
# goodm, goodmDown: whether the Goodman et al. 1995 restriction barring multiple-dose jumps is in effect, 
#               for the up ('goodm') and down ('goodmDown') directions.
# ...: passed on to dfcrm::crm(). In particular, the scale argument controls the prior distribution's SD.
{
require(dfcrm)
tmpp = crm(prior=skel, target=targ, tox=responses, level=doses, var.est=FALSE, patient.detail=FALSE, model.detail=FALSE)
currdose = doses[length(doses)]

# Goodman et al. escalation constraint
if (goodm) cand=min(tmpp$mtd,currdose+1)  
# Also in the down direction...
if (goodmDown && tmpp$mtd<currdose) cand=currdose-1   

if (pace>1 && cand>currdose)  ## install safety pace for each escalation
{
	runs=rle(doses)
	if(runs$lengths[length(runs$lengths)] < pace) cand=currdose
}
return(cand)
}

#### Ivanova et al.'s Cumulative Cohort Design (CCD), one of the first interval designs

### wrapping CCD

ccd <- function(doses, responses, hwidth=0.1, targ=0.3, pace=1, shrink=FALSE, swt=1, boundout = c(TRUE, TRUE), ...)
{
n=length(doses)
if(length(responses)!=n) stop("Doses, reponses vector unequal in length.")
currdose = doses[n]

testy=mean(responses[doses==currdose])
if(shrink) {
	nm=sum(doses==currdose)
	testy=(testy+swt*targ)/(nm+swt)
}

# Above interval; go down
if(testy > targ+hwidth) return(currdose-1)
# Boundary case
if(testy == targ+hwidth && boundout[1]) return(currdose-1)
# Below interval: go up
if(testy < targ-hwidth || (testy == targ-hwidth && boundout[2]) )
{
	if (pace>1)  ## optional safety pace for each escalation
	{
		runs=rle(doses)
		if(runs$lengths[length(runs$lengths)] < pace) return(currdose)
	}
	return(currdose+1)
}
# Otherwise...... stay:
return(currdose)
}

#### Implementing Liu and Yuan (2015)'s BOIN interval design.
# Package not required for this function, but look-up table such as that produced via BOIN::get.boundary() is needed

boin <- function(doses, responses, lookup, hardstop = FALSE, pace = 1, excmin = 3, ...)
{
n=length(doses)
if(length(responses)!=n) stop("Doses, reponses vector unequal in length.")
if(any(diff(lookup[1,]) > 1) ) stop("We prefer the lookup table to have increments of 1.")
if(max(lookup[1,]) < n)  stop("Lookup table too short.")
currdose = doses[n]
# Getting current dose's n and y
currn = sum(doses == currdose)
curry = sum(responses[doses == currdose])

### The lookup table's first row is the current n (i.e., the lookup value)
refnums = lookup[2:4, lookup[1,] == currn] 

# De-escalation, plus checking whether entire trial is shut down
if(curry >= refnums[2]) {
	if(curry >= refnums[3] && currdose==1 && hardstop) return(NA)
	return(currdose - 1)
}

# Before escalating, we must check next dose still "open for business"
if(curry <= refnums[1]) 
{
	if (pace>1)  ## optional safety pace for each escalation
	{
		runs=rle(doses)
		if(runs$lengths[length(runs$lengths)] < pace) return(currdose)
	}
	nextn = sum(doses == currdose + 1)
	if(nextn < excmin) return(currdose + 1)
	nextexc = lookup[4, lookup[1,] == nextn] 
	if(is.na(nextexc)) return(currdose + 1)
	nexty = sum(responses[doses == currdose + 1])
	if(nexty < nextexc) return(currdose + 1)
}
return(currdose)
}

#### Lastly... the world-famous 3+3 !
# Note that this is only the escalation/stopping protocol. 
# MTD selection requires a separate function.

threePlus3 <- function(doses, responses, maxout, minout=1, noReEsc=TRUE, direc=FALSE,...)
{
n=length(doses)
curr=doses[n]
if(direc) curr=0
if(n%%3 > 0) return(curr) # only evaluating after threes
if(n<3) return(curr)

currn=which(doses==doses[n])
if (length(currn)==3) {
	dlt=sum(responses[(n-2):n])
	if (dlt==0 && noReEsc && doses[n]<max(doses)) return(curr)  ### repeating due to never-re-escalate rule
	dout=curr+sign(1-dlt)
} else {
	dlt=sum(responses[currn])
	if (dlt<2 && noReEsc && doses[n]<max(doses)) return(NA)  ### stopping due to never-re-escalate rule
	dout=curr+sign(1.5-dlt)
}
if (ifelse(direc,doses[n]+dout,dout)>maxout) dout=maxout
if (ifelse(direc,doses[n]+dout,dout)<minout) dout=minout

### Stopping rules?
outn=length(doses[doses==ifelse(direc,doses[n]+dout,dout)])
if(outn==6) return(NA)
return(dout)
}


