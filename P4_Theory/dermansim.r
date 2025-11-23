rm(list=ls())
library(magrittr)
library(data.table)

pweib3 <- function(x, shp, scl, shift) pweibull(q=x+shift, shape=shp, scale=scl)
qweib3 <- function(p, shp, scl, shift) qweibull(p, shape=shp, scale=scl) - shift



derman <- function(doses, responses, coin, lowTarget, fastStart=FALSE,...)
{
  n=length(doses)
  curr=doses[n]
  if(!lowTarget) {
    dout=ifelse(responses[n]==0,curr+1,ifelse(runif(1)<=coin,curr-1,curr+1))
    if(fastStart && sum(responses)==n) dout=curr-1
  } else {
    dout=ifelse(responses[n]==1,curr-1,ifelse(runif(1)<=coin,curr+1,curr-1))
    if(fastStart && sum(responses)==0) dout=curr+1
  }
  return(dout)
}


outdir = 'C:/GitHub/output'
weib30parm = fread(file.path(outdir, 'scenarios_weib30.csv') )
weib90parm = fread(file.path(outdir, 'scenarios_weib90.csv') )

set.seed(179286)
nsim = 500 # This overrides 'nsim' in the scenario-generation scripts
M = 8
n = 50
b30list = list(coin=3/7, lowTarget=TRUE, fastStart=FALSE)
d30list = list(coin=5/7, lowTarget=TRUE, fastStart=FALSE)
b90list = list(coin=1/9, lowTarget=FALSE, fastStart=TRUE)
d90list = list(coin=5/9, lowTarget=FALSE, fastStart=TRUE)

lostart = 1

# Generate the scenarios on the dose grid
weib30F = weib30parm[1:nsim, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
					
# Generate response thresholds ("patients")
#  (we use the same thresholds throughout, to eliminate this variability source)

thresh30w = matrix(runif(n*nsim), nrow=n)

bw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						design=bcd, desArgs=b30list, thresholds = thresh30w)					
dw30minmid = dfsim(n, starting = lostart, Fvals=weib30F, ensemble=nsim,
						design=derman, desArgs=d30list, thresholds = thresh30w)					

M9 = 12
weib90F = weib90parm[1:nsim, ][ ,as.vector(pweib3(1:M9, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M9)

bw90minmid = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh30w)					
dw90minmid = dfsim(n, starting = lostart, Fvals=weib90F, ensemble=nsim,
						design=derman, desArgs=d90list, thresholds = thresh30w)					

bw90maxmid = dfsim(n, starting = M9, Fvals=weib90F, ensemble=nsim,
						design=bcd, desArgs=b90list, thresholds = thresh30w)					
dw90maxmid = dfsim(n, starting = M9, Fvals=weib90F, ensemble=nsim,
						design=derman, desArgs=d90list, thresholds = thresh30w)					
					