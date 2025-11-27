### 
#------------------ Utility: coin calculation for randomized GUD
#  Assaf Oron, Fall 2025

### The utility provides the coin probability to produce a desired target.
### The probability returned is for the "higher" of the two options.

#### Arguments:
#  cohort, lower, upper: GUD parameters, similar to the upndown function gudtarg()
#  randomize: the outcome to randomize on. At present must equal the value of either lower or upper.
#  target: the desired exact target response rate between 0 and 1.
#  digits: the rounding resolution.

#### NOTE! If the design combination you specified cannot yield the desired target, you will receive
#      the standard R uniroot() error message: "f() values at end points not of opposite sign".
#      In that case, check again by running gudtarg() on the two non-randomized cases.
#      For example, if using cohort = 3, lower = 0, upper = 1, and randomize = 1
#       check gudtarg(3,0,1) and gudtarg(3,0,2), and see that your desired target falls between them.

gudcoin <- function(cohort, lower, upper, randomize, target, digits = 4)
{
require(upndown)
checkTarget(target)
checkNatural(c(cohort, lower+1, upper), parname = 'cohort, lower+1, upper', toolarge = 50)  
if(cohort<upper || upper<=lower) stop('Order must be lower < upper <= cohort.\n')
if(!(randomize %in% c(lower, upper)) ) stop('The value to be randomized must be on the decision boundary.\n')

# Cases
if(randomize == upper)
{
	b = uniroot(f=function(x, tar, k, u, l) {pbinom(q=l, size=k, prob=tar) + 
			ifelse(u==k, 1, pbinom(q=u, size=k, prob=tar) ) - 1 - (1-x) * dbinom(x=u, size=k, prob=tar)}, 
			interval=0:1, k=cohort, u=upper, l=lower, tar=target)$root
}  else if(randomize == lower)
{
	b = uniroot(f=function(x, tar, k, u, l) {ifelse(l==0, 0, pbinom(q=l-1, size=k, prob=tar) ) + 
			pbinom(q=u-1, size=k, prob=tar) - 1 + x * dbinom(x=l, size=k, prob=tar)}, 
			interval=0:1, k=cohort, u=upper, l=lower, tar=target)$root
} 
cat("After Y = ", randomize, ", with probability", round(b, digits), ifelse(randomize==lower, "go UP", "REPEAT"), '\n')
cat("          and with probability", round(1-b, digits), ifelse(randomize==lower, "REPEAT.", "go DOWN."), '\n')

}