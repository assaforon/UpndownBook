liao_PAVA = function(y, n) #crashes when any element of n is 0
{
  ratio = cumsum(y) / cumsum(n)
  k_min = ratio |> which.min()
  
  isotonic_est = rep(ratio[k_min], k_min)
  
  m = length(y)
  if(k_min == m) isotonic_est
  else 
  {
    index1 = (k_min+1):m
    c(isotonic_est, liao_PAVA(y[index1], n[index1])) 
  }
}
########
log_likelihood_partial = function(yy, nn, theta_star) #only looks up
{   
  #H_0: theta_1 = theta_star vs. H_1: \theta_1 < theta_star
  mm = yy |> length()
  ##############################
  theta_free = liao_PAVA(yy, nn)  
  if(theta_free[1] >= theta_star) return(0)
  
  index1 = theta_free < theta_star
  
  theta_free = theta_free[index1]
  yy = yy[index1]
  nn = nn[index1]
  
  small = 1.e-6
  theta_free[theta_free < small] = small
  theta_free[theta_free > 1 - small] = 1 - small
 
  ################
  value1 = yy*log(theta_free/theta_star) + (nn - yy)*log( (1-theta_free)/(1-theta_star))
  value1 = sum(value1)
  
  2*sum(value1) # log-likelihood, 
}

##################################
my_uniroot = function(func1, lower, upper) #func1 is decreasing
{
  f.lower = func1(lower)
  if(f.lower <= 0) return(lower)
  f.upper = func1(upper)
  if(f.upper >= 0) return(upper)
  
  uniroot(func1, lower=lower, upper=upper, f.lower=f.lower, f.upper=f.upper, maxiter = 30)$root
}

rbinom_quasi = function(size, prob, n_simu)
{
  q = seq(0, 1, length=n_simu+2)
  qbinom(q[2:(n_simu+1)], size, prob) |> sample()
}

p_value_less_partial = function(y, n, i, theta_star, s) #alternative is theta_i < theta_star
{ 
  m = length(y)
  upper = min(i + s - 1, m)
  yy = y[i:upper]
  nn = n[i:upper]
  
  observed = log_likelihood_partial(yy, nn, theta_star)
  
  n_simu = 2000
  draws = vapply(nn, FUN=rbinom_quasi, FUN.VALUE=numeric(n_simu), n_simu = n_simu, prob=theta_star)
  #draws n_simu by length(nn)
  simulated = apply(draws, 1, FUN=log_likelihood_partial, nn=nn, theta_star=theta_star)
  
  mean(simulated > observed) + mean(simulated == observed)/2 #mid-pvalue
}


#######################################
one_upper_limit_partial = function(y, n, i, conf.level_one_sided, s)
{ 
  f1 = function(theta_star) p_value_less_partial(y, n, i, theta_star, s) - (1 - conf.level_one_sided)
  
  my_uniroot(f1, lower=0.0001, upper=0.9999)
}

one_lower_limit_partial = function(y, n, i, conf.level_one_sided, s)
{ 
  m = length(y)
  1 - one_upper_limit_partial(rev(n - y), rev(n), m-i+1, conf.level_one_sided, s)
}


###########################################

liaoCI = function(y, n, phat = y/n, narrower = TRUE, alternate=wilsonCI,
				conf = 0.9, s = 3, ...) #conf.level is the two-sided confidence
  #s is the number of dose-levels to try to combine information 
{
  #' @param y integer or numeric vector, the pointwise Binomial counts
  #' @param n integer or numeric vector, the pointwise sample sizes
    
  if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
  
  conf.level_one_sided = 1 - (1 - conf)/2
  
  f1 = function(i)
  {
    term1 = one_lower_limit_partial(y, n, i, conf.level_one_sided, s)
    term2 = one_upper_limit_partial(y, n, i, conf.level_one_sided, s)
    c(term1, term2)
  }
  
  m = length(y)
  cand = sapply(1:m, f1) |> t()
  
  lcl = cand[ ,1]
  ucl = cand[ ,2]
  
   if(narrower) # Optional pointwise narrowing via alternate interval function
# See Oron and Flournoy (2017) for justification/performance
{
	relevants=which(n>0) # Avoiding n=0 boundary where alternate() produces NaNs
	altout=alternate(phat=phat[relevants], n=n[relevants], conf=conf, ...)
# The cummax, cummin (added Dec. 2015) ensure monotonicity of boundaries.
# Monotonicity is not for the optics, but rather another way to pool information
# from where it is plentiful to where it might be lacking.

	cand[relevants, 1] = cummax(pmax(lcl[relevants],altout[,1]))
	cand[relevants, 2] = rev(cummin(rev(pmin(ucl[relevants],altout[,2]))))
}   
 cand
}



two_one_side_p_values = function(y, n, delta1, delta2)
{
  m = y |> length()
  
  f1 = function(i) 
  {
    term1 = p_value_greater_partial(y[1:i], n[1:i], delta1)
    term2 = p_value_less_partial(y[i:m], n[i:m], delta2)
    c(term1, term2)
  }
  
  cat("p-value for resting H_0: theta <= ", delta1, "and H_0: theta >= ", delta2, "\n")
  p_values = sapply(1:m, f1) |> t()
  p_values |> print()
  p_values
}

