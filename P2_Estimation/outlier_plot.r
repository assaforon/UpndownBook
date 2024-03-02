cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

load(file.path(outdir, 'grandsim50l.RData'))	

pdf(file.path(outdir, 'outlier_plot.pdf'), width = 11, height = 6)
xtit = 'Simulation Run (reordered by error size)'
ytit = 'Square Estimation Error'
ymax = estl50lomid[ , max((cir-true)^2, na.rm = TRUE) ]

layout(t(1:2), widths = 20:19)
par(mar = c(3, 3., 3, 1), mgp = c(1.7, 0.4, 0), las = 1, xaxt = 'n', tck = -0.01)
psize = 0.5

estl50lomid[ , plot(sort((cir-true)^2), xlab = xtit, ylab = ytit, cex = psize,
	main = 'CIR', ylim = c(0, ymax) ) ]
abline(h = 1, lty = 2)
	
	par(mar = c(3, 1, 3, 1) )
estl50lomid[ , plot(sort((dyna-true)^2), xlab = xtit, ylab = '', cex = psize,
	main = 'Dynamic Average', ylim = c(0, ymax) ) ]
abline(h = 1, lty = 2)

dev.off()
cat(base::date(), '\n')
