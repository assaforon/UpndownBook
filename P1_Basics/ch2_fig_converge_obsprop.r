source('basics_header.r')

longdat = readRDS(file.path(outdir, 'ch2_sim_converge_obsprop.rds') )
longn = nrow(longdat$doses)
nsim = dim(longdat$doses)[2]

### Additional constants/utilities
bpi = pivec(exampleF, bcdmat, target = btarg)

# Cumulative allocation frequency of dose m
cumulprop = function(x, m) cumsum(x==m)/(1:length(x)) 
runningprop = function(x, tstart, slice, m)  mean(x[tstart:(tstart + slice -1)] == m)


p4 = apply(longdat$doses, 2, cumulprop, m=4)
p410 = apply(longdat$doses[11:longn, ], 2, cumulprop, m=4)

#run4 = matrix(NA, nrow = longn - 20, ncol = nsim)
#for(a in 1:(1200) ) run4[a, ] = apply(longdat$doses, 2, runningprop, tstart = a, slice = 200, m=4)
p7 = apply(longdat$doses, 2, cumulprop, m=7)

pivecs = matrix(NA, nrow = 100, ncol = M)
for(a in 1:100) pivecs[a,] = currentvec(exampleF, bcdmat, n=a, startdose=1, target = 0.3)

plotstart = 10
plotend = 1100
xlog = log(plotstart:plotend)

layout(t(1:2))
par(tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)

plot(log(abs(rowMeans(p4[plotstart:plotend,])-bpi[4])) ~ xlog, cex = 0.5, xaxt = 'n', yaxt = 'n', xlim = log(c(plotstart, 1000) ),
	ylim = log(c(10^-4, 0.15)),
	main = expression(paste('BCD Example, ', d[4]) ), xlab = 't', ylab = expression(paste('SD, or Absolute Difference from  ', pi) ) )
xvals = sort(as.vector(outer(c(1,3), 10^(1:3) ) ) )
axis(1, at = log(xvals), labels =  xvals)
yvals = sort(as.vector(outer(c(1,3), 10^-(1:4) ) ) )
axis(2, at = log(yvals), labels =  yvals )

points(log(apply(p4[plotstart:plotend,], 1, sd)) ~ xlog, pch = 'I', cex = 0.5)
points(log(abs(rowMeans(p410[1:(plotend-9),])-bpi[4])) ~ log(plotstart:plotend), pch='x', cex = 0.5 )
lines(log(abs(pivecs[,4]-bpi[4])) ~ log(1:100) )


plot(log(abs(rowMeans(p7[plotstart:plotend,])-bpi[7])) ~ xlog, cex = 0.5, xaxt = 'n', yaxt = 'n', xlim = c(log(c(plotstart, 1000) ) ),
	main = expression(paste('BCD Example, ', d[7]) ), xlab = 't', ylab = expression(paste('SD, or Absolute Difference from  ', pi) ) )
axis(1, at = log(xvals), labels =  xvals)
axis(2, at = log(yvals), labels =  yvals )

points(log(apply(p7[plotstart:plotend,], 1, sd)) ~ xlog, pch = 'I', cex = 0.5)
lines(log(abs(pivecs[,7]-bpi[7])) ~ log(1:100) )



