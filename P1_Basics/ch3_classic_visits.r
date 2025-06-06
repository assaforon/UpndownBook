source('basics_header.r')

startdose = 2

cpi = pivec(exampleF, classicmat)

eye = diag(1, M, M) # identity matrix I
Z = solve(eye - classicmat(exampleF) + t(matrix(rep(cpi, M), nrow = M)) )
Zdg = diag(diag(Z)) # version of Z with 0 off-diags
# K & S book matrices
dee = diag(1/cpi, M, M) # D
E = matrix(1, M, M) 

initvec = rep(0, M)
initvec[startdose] = 1

# Mean first-passage time per Kemeny and Snell
Mij =  (eye - Z + E %*% Zdg ) %*% dee
fij = as.vector(initvec %*% Mij)

# Second moment and SD (not used in plot):
Wij = Mij %*% (2 * Zdg %*% dee - eye) + 2* (Z %*% Mij - E %*% diag(diag(Z %*% Mij)) )
f2ij = as.vector(initvec %*% Wij)

fsd = sqrt(f2ij - fij^2)

######## Simulated distribution of visits

classic5k = readRDS(file.path(outdir, 'ch3_sim_firstvisit_distribution_classical.rds') ) 
longn = nrow(classic5k$doses)
nsim = ncol(classic5k$doses)

# Calculating a matrix of first visits (m in rows, run ID in columns)
firstvisit = function(x, m) min(which(x == m))

firsts = matrix(NA, nrow = M, ncol = nsim )
for(m in 1:M) firsts[m,] = apply(classic5k$doses[-1, ], 2, firstvisit, m=m)

# 90% simulation interval
visit90 = apply(firsts, 1, quantile, prob = c(1,19)/20, type = 6)


########## Plotting

yseq = seq(0,50,10)

pdf(file.path(outdir, 'ch3_classic_visits.pdf'), width = 17, height = 6)

layout(t(1:3))
par(tck = -0.01, las = 1, cex.lab = 2, cex.axis = 1.7, cex.main = 2.5, mar = c(4.,4.5,4,1), mgp = c(3,0.5,0) )

# Distribution
plot(c(1,M), log10(c(1,1500)), type='n', xaxt='n', yaxt='n', main = expression(paste('First Visit from ', 
					d[2], ': Mean and Uncertainty') ), xlab = 'Destination Dose-Level', ylab = 'Observations until First Visit')

polygon(x=c(1:M, M:1), y=log10(c(visit90[1,], rev(visit90[2,]))), col='grey85', border='grey85')
lines(1:M, log10(fij), lwd = 2)

axis(1, 1:M)
ynums = sort(outer(c(1,3), 10^(0:3)) )
axis(2, at = log10(ynums), labels = ynums)

abline(v = 5.6, lty = 2)

# Demo of distribution from d_2 to d_6, well-behaved

plot(table(firsts[6,]),xaxt='n', yaxt='n', main = expression(paste('First Visit from ', 
					d[2], ' to ', d[6], ': Full Distribution' ) ), xlab = 'Observations until First Visit', 
					ylab = paste('Percent of Simulated Experiments'), xlim = c(4,35))
axis(1, at=seq(5,40,5) )
axis(2, at=yseq, labels = yseq/2)
abline(v=fij[6], lty = 2)
abline(v=mean(firsts[6,]), lty = 3)

legend('topright', legend = paste(c('Theoretical', 'Simulation'), 'Mean'), lty = 2:3, bty = 'n', cex = 1.8)

# Demo of distribution from d_2 to d_1, particularly skewed
plot(table(firsts[1,]),xaxt='n', yaxt='n', main = expression(paste('First Visit from ', 
					d[2], ' to ', d[1], ': Full Distribution' ) ), xlab = 'Observations until First Visit', 
					ylab = paste('Percent of Simulated Experiments'), ylim = c(0, 49) )
axis(1, at=seq(0,3000,500) )
axis(2, at=yseq, labels = yseq/2)
abline(v=fij[1], lty = 2)
abline(v=mean(firsts[1,]), lty = 3)




dev.off()
