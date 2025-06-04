source('basics_header.r')

### Additional constants/utilities
wid = 2

# Stationary distributions
bpi90 = 100*pivec(exampleF2, bcdmat, target = 0.9)
bpi95 = 100*pivec(exampleF2, bcdmat, target = 0.95)

showrange = 5:17

#-------------  Simple F(x), pi figure of the two distributions
pdf(file.path(outdir, 'ch3_bcd_fig9095.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

p1 = barplot(bpi90[showrange], 
	xlab = 'Dose-Level', main = 'Targeting the 90th Percentile', ylim = c(0, max(bpi95)),
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ) )
axis(1, at = p1, labels = showrange)

p2 = barplot(bpi95[showrange], xlab = 'Dose-Level', main = 'Targeting the 95th Percentile',
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ) )
axis(1, at = p2, labels = showrange)

# abline(v = qlogis(0.9, location = exampleMu, scale = exampleSig), lty=2, lwd=wid) 

# abline(v = qlogis(0.95, location = exampleMu, scale = exampleSig), lty=2, lwd=wid) 

dev.off()

