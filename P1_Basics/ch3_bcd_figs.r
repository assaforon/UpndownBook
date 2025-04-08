source('basics_header.r')

### Additional constants/utilities

startdose = 2
initvec = rep(0, M)
initvec[startdose] = 1
btarg = 0.3

# Stationary distributions
bpi = pivec(exampleF, bcdmat, target = btarg)
ebpi = pivec(exp11F, bcdmat, target = btarg)

#-------------  Simple F(x), pi figure of the two distributions
pdf(file.path(outdir, 'ch3_bcd_pi.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

plot(exampleF, type = 'l', lwd=2, xaxt = 'n', xlab = 'Dose-Level', ylab = 'F(x)', ylim = c(0.03,0.9) )
lines(exp11F, lty=2, lwd=2)
abline(h = 0.3, lty = 3)
axis(1, at = 1:M)

plot(100 * bpi, type = 'l', lwd=2, xaxt = 'n', xlab = 'Dose-Level', 
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), ylim = c(1, 27) )
lines(100 * ebpi, lty=2, lwd=2)
axis(1, at = 1:M)
legend('topright', legend = c('Logistic', 'Exponential'), lty = 1:2, lwd = 2, cex = 1.2, bty = 'n')

dev.off()

#------------ Visit times

# Matrices for visit calculation

eye = diag(1, M, M) # identity matrix I
E = matrix(1, M, M) 

# Specific matrices: logistic curve
Zl = solve(eye - bcdmat(exampleF, target = btarg) + t(matrix(rep(bpi, M), nrow = M)) )
Zldg = diag(diag(Zl)) # version of Z with 0 off-diags
deel = diag(1/bpi, M, M) # D

# Specific matrices: exponential curve
Ze = solve(eye - bcdmat(exp11F, target = btarg) + t(matrix(rep(ebpi, M), nrow = M)) )
Zedg = diag(diag(Ze)) # version of Z with 0 off-diags
deee = diag(1/ebpi, M, M) # D

### Mean first-passage time per Kemeny and Snell
# logistic
Mlij =  (eye - Zl + E %*% Zldg ) %*% deel
flij = as.vector(initvec %*% Mlij)
#exponential
Meij =  (eye - Ze + E %*% Zedg ) %*% deee
feij = as.vector(initvec %*% Meij)

### Simulated distributions for variability bands

longdat1 = readRDS(file.path(outdir, 'ch3_sim_firstvisit_distribution_bcd1.rds') )
longdat2 = readRDS(file.path(outdir, 'ch3_sim_firstvisit_distribution_bcd2.rds') )
longn = nrow(longdat1$doses)
nsim = ncol(longdat1$doses)

# Calculating a matrix of first visits (m in rows, run ID in columns)
firstvisit = function(x, m) min(which(x == m))

lfirsts = matrix(NA, nrow = M, ncol = nsim )
efirsts = lfirsts
for(m in 1:M) 
{
	lfirsts[m,] = apply(longdat1$doses[-1, ], 2, firstvisit, m=m)
	efirsts[m,] = apply(longdat2$doses[-1, ], 2, firstvisit, m=m)
}

# 90% simulation interval
lvisit90 = apply(lfirsts, 1, quantile, prob = c(1,19)/20, type = 6)
evisit90 = apply(efirsts, 1, quantile, prob = c(1,19)/20, type = 6)

### Plotting constants 

bands = 'grey85'
ynums = sort(outer(c(1,3), 10^(0:4)) )
yrange = c(1, 1e4)
lvisit90[!is.finite(lvisit90)] = 2 * yrange[2]

pdf(file.path(outdir, 'ch3_bcd_visits.pdf'), width = 14, height = 7.3)

layout(t(1:2), widths = 16:15)
par(tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7, mar = c(4.5,4.5,4,1), mgp = c(3,0.5,0) )

# logistic
plot(c(1,M), log10(yrange), type='n', xaxt='n', yaxt='n', main = expression(paste('First Visit from ', 
					d[2], ', Logistic F(x)') ), xlab = 'Destination Dose-Level', ylab = 'Observations until First Visit')

polygon(x=c(1:M, M:1), y=log10(c(lvisit90[1,], rev(lvisit90[2,]))), col=bands, border=bands)
lines(1:M, log10(flij), lwd = 2)

axis(1, 1:M)
axis(2, at = log10(ynums), labels = ynums)
abline(v = 3.9, lty = 2)

# exponential
par(mar = c(4.5,3,4,1))
plot(c(1,M), log10(yrange), type='n', xaxt='n', yaxt='n', main = expression(paste('First Visit from ', 
					d[2], ', Exponential F(x)') ), xlab = 'Destination Dose-Level', ylab = '')

polygon(x=c(1:M, M:1), y=log10(c(evisit90[1,], rev(evisit90[2,]))), col=bands, border=bands)
lines(1:M, log10(feij), lwd = 2)

axis(1, 1:M)
axis(2, at = log10(ynums), labels = ynums)
abline(v = 3.9, lty = 2)

dev.off()



