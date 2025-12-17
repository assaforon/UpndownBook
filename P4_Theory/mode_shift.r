source('../P1_Basics/basics_header.r')

### Additional constants/utilities
lwid = 1.5
xshift = (1:M) + 0.5
btarg = 0.3

exampleFshift = plogis(xshift, location = exampleMu, scale = exampleSig)


# Stationary distributions
c0 = 100*pivec(exampleF, classicmat)
c5 = 100*pivec(exampleFshift, classicmat)

b0 = 100*pivec(exampleF, bcdmat, target = btarg)
b5 = 100*pivec(exampleFshift, bcdmat, target = btarg)

showrange = 5:17

#-------------  Simple F(x), pi figure of the two distributions
pdf(file.path(outdir, 'mode_shift.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

plot(c0, type = 'l', lwd = lwid, xlim = c(1, M+0.5), ylim = c(1, max(c(c0,c5)) ), main = 'Classical UDD',
	xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)', xaxt = 'n' ) 
lines(xshift, c5, lwd = lwid, lty = 2)
points(xshift, c5)
points(c0, pch = 19)
axis(1, at = 1:M)
abline(v = exampleMu, lty = 3)

plot(b0, type = 'l', lwd = lwid, xlim = c(1, M+0.5), ylim = c(1, max(c(b0,b5)) ), main = 'BCD (30th percentile)', 
	xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)', xaxt = 'n' ) 
lines(xshift, b5, lwd = lwid, lty = 2)
points(xshift, b5)
points(b0, pch = 19)
axis(1, at = 1:M)
abline(v = qlogis(btarg, location = exampleMu, scale = exampleSig), lty = 3)

dev.off()


