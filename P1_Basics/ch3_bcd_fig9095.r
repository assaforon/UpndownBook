source('basics_header.r')

### Additional constants/utilities
wid = 2

# Stationary distributions
bpi90 = 100*pivec(exampleF2, bcdmat, target = 0.9)
bpi95 = 100*pivec(exampleF2, bcdmat, target = 0.95)

#-------------  Simple F(x), pi figure of the two distributions
pdf(file.path(outdir, 'ch3_bcd_fig9095.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)


plot(bpi90, type = 'l', lwd=wid, xaxt = 'n', xlab = 'Dose-Level', main = 'Targeting the 90th Percentile',
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), ylim = c(1, max(bpi95)), xlim=c(5,17) )
abline(v = qlogis(0.9, location = exampleMu, scale = exampleSig), lty=2, lwd=wid) 
axis(1, at = 1:(2*M) )

plot(bpi95, type = 'l', lwd=wid, xaxt = 'n', xlab = 'Dose-Level', main = 'Targeting the 95th Percentile',
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), ylim = c(1, max(bpi95)), xlim=c(5,17)  )
abline(v = qlogis(0.95, location = exampleMu, scale = exampleSig), lty=2, lwd=wid) 
axis(1, at = 1:(2*M) )

dev.off()

