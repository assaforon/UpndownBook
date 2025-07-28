source('basics_header.r')

longdat = readRDS(file.path(outdir, 'ch2_sim_converge_obsprop.rds') )
longn = nrow(longdat$doses)
nsim = dim(longdat$doses)[2]

### Additional constants/utilities
bpi = pivec(exampleF, bcdmat, target = btarg)

# Cumulative allocation frequency of dose m
cumulprop = function(x, m) cumsum(x==m)/(1:length(x)) 
runningprop = function(x, tstart, slice, m)  mean(x[tstart:(tstart + slice -1)] == m)

# Cumulative at d_4, closest to target
p4 = apply(longdat$doses, 2, cumulprop, m=4)
# Cumulative at d_7, ~11x less often and far from starting point
p7 = apply(longdat$doses, 2, cumulprop, m=7)

pivecs = matrix(NA, nrow = 999, ncol = M)
for(a in 1:999) pivecs[a,] = currentvec(exampleF, bcdmat, n=a, startdose=1, target = 0.3)

### Plotting parametersj
plotstart = 11
plotend = 1100
psize = 0.9
lwid = 2

xlog = log(plotstart:plotend)

pdf(file.path(outdir, 'ch2_converge_obsprop.pdf'), width = 14, height = 7.3)

layout(t(1:2), widths = 16:15 )
par(tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7, mar = c(5,6,4,1), mgp = c(3.5,0.5,0) )

plot(log(abs(rowMeans(p4[plotstart:plotend,])-bpi[4])) ~ xlog, cex = psize, xaxt = 'n', yaxt = 'n', xlim = log(c(plotstart+0.5, 999) ),
	ylim = log(c(5e-4, 0.15)),
	main = expression(paste('BCD Example, Allocations to ', d[4]) ), xlab = 't (cumulative observations)', ylab = expression(paste('SD, or Absolute Difference from  ', pi) ) )
xvals = sort(as.vector(outer(c(1,3), 10^(1:3) ) ) )
axis(1, at = log(xvals), labels =  xvals)
yvals = sort(as.vector(outer(c(1,3), 10^-(1:4) ) ) )
axis(2, at = log(yvals), labels =  yvals )

points(log(apply(p4[plotstart:plotend,], 1, sd)) ~ xlog, pch = 'I', cex = psize)
lines(log(abs(cumsum(pivecs[,4])/(1:999)-bpi[4])) ~ log(1:999), lwd = lwid )
lines(log(abs(pivecs[,4]-bpi[4])) ~ log(1:999), lty = 2 )

legend('bottom', legend = c('p(t)', 'Cumulative p(t)', paste('N(t)/t:', c('Mean', 'SD')) ), pch = c(NA, NA, 'o', 'I'), 
				lwd = c(1, lwid, 0, 0), lty = c(2,1,0,0), bty = 'n', cex = 1.3, pt.cex = 1.2)


par(mar = c(5,3,4,1) )
plot(log(abs(rowMeans(p7[plotstart:plotend,])-bpi[7])) ~ xlog, cex = psize, xaxt = 'n', yaxt = 'n', 
	xlim = c(log(c(plotstart+0.5, 999) ) ), ylim = log(c(3e-4, 0.035)),
	main = expression(paste('BCD Example, Allocations to ', d[7]) ), xlab = 't (cumulative observations)', ylab = '' )
axis(1, at = log(xvals), labels =  xvals)
axis(2, at = log(yvals), labels =  yvals )

points(log(apply(p7[plotstart:plotend,], 1, sd)) ~ xlog, pch = 'I', cex = psize)
lines(log(abs(cumsum(pivecs[,7])/(1:999)-bpi[7])) ~ log(1:999), lwd = lwid )
lines(log(abs(pivecs[,7]-bpi[7])) ~ log(1:999), lty = 2 )


dev.off()


