source('basics_header.r')

n = 40
nsim = 3
M = length(exampleF)

set.seed(2025) # Comment out this line, to get a different ensemble each time

thresh = matrix( runif(n*nsim), nrow = n)

# Simulating the runs

classic = dfsim(n = n, starting = M/2, Fvals = exampleF, ensemble = nsim, 
	design = krow, desArgs = list(k=1), thresholds = thresh )
p80 = dfsim(n = n, starting = M/2, Fvals = exampleF, ensemble = nsim, 
	design = krow, desArgs = list(k=3, lowTarget = FALSE), thresholds = thresh )
p30 = dfsim(n = n, starting = M/2, Fvals = exampleF, ensemble = nsim, 
	design = bcd, desArgs = list(k=3, coin = 3/7, lowTarget = TRUE), thresholds = thresh )
	
# Targets for showing on the plots
t50 = qlogis(0.5, location = exampleMu, scale = exampleSig)
t80 = qlogis(k2targ(3, lowTarget = FALSE), location = exampleMu, scale = exampleSig)
t30 = qlogis(0.3, location = exampleMu, scale = exampleSig)

# Plotting

hstyle = 3
psize = 1.4
shp = 'square'

pdf(file.path(outdir, 'ch1_randomruns.pdf'), width = 15, height = 13)
layout(matrix(1:9, nrow = 3), heights = c(18,18,19), widths = c(15,14,14) )
par(mgp = c(2.2, 0.6, 0), tck = -0.01, las = 1, cex.lab = 2, cex.axis = 1.4, cex.main = 2)

for(a in 1:nsim)
{
	lmar = ifelse(a>1, 1., 4)
	par(mar = c(2.5,lmar,4.,1.) )
	yl = ifelse(a==1, 'Dose-Level', '')
	udplot(classic$doses[ ,a], classic$responses[ ,a], allow1extra = TRUE, ylim = c(1, M), 
		cex = psize, yaxt = 'n', xtitle = '', ytitle = yl, shape = shp, main = 'Classical' )
	axis(2, 1:M); abline(h = t50, lty = hstyle, lwd = 1.5)
	udplot(p80$doses[ ,a], p80$responses[ ,a], allow1extra = TRUE, ylim = c(1, M), 
		cex = psize, yaxt = 'n', xtitle = '', ytitle = yl, shape = shp, main = 'K-in-a-row' )
	axis(2, 1:M); abline(h = t80, lty = hstyle, lwd = 1.5)
	par(mar = c(4,lmar,4.,1.) )
	udplot(p30$doses[ ,a], p30$responses[ ,a], allow1extra = TRUE, ylim = c(1, M), 
		cex = psize, yaxt = 'n', ytitle = yl, shape = shp, main = 'Biased-Coin')
	axis(2, 1:M); abline(h = t30, lty = hstyle, lwd = 1.5)
}

dev.off()

