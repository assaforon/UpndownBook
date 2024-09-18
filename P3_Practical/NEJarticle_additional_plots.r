library(upndown)

#--------------------- CIR plot

dlabel = 'Phenylephrine dose (micrograms)'
george10x = 80 + 20 * c(1, rep(2, 5), 1, 1, 0, 0, rep(1, 7), 0:2, 2, 2, rep(1, 4), 2, 1, 1, 2, 2, 
                        rep(3, 5), 4, 5, 5, rep(4, 6))
george10y = c(ifelse(diff(george10x) > 0, 0, 1), 1)

pdf('../../output/NEJarticle_CIR.pdf', width = 9, height = 8)

par(las = 1, cex.axis = 1.3, cex.lab = 1.6, mar = c(4.5,4.5,1,1) )
drplot(x=george10x, y=george10y, addest = TRUE, target = 0.9, addcurve = TRUE, 
			balancePt = 10/11, xtitle = dlabel, percents = TRUE, ytitle = "Efficacy (%)" )

irest = quickIsotone(x=george10x, y=george10y, estfun = oldPAVA)

lines((100*y)~x, data=irest, lty=2)

legend('bottomright', bty = 'n', lty = c(0,2,1,1), pch =  c(4,NA,NA,19), col = c(1,1,'blue','purple'),
				legend = c("Observations", 'IR Curve', 'CIR Curve', 'Target Estimate'), cex = 1.4, lwd = 2)
				
dev.off()

