source('basics_header.r')


btarg = 0.3

cpi = pivec(exampleF, classicmat)
bpi = pivec(exampleF, bcdmat, target = btarg)


pdf(file.path(outdir, 'ch2_pi.pdf'), width = 14, height = 7.3)

layout(t(1:2))

par(mar = c(4,4,4,2), mgp = c(2.5, 0.6, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)
# We need to assign the plot, b/c barplot() returns the true x values at the mid-bars
p1 = barplot(100*cpi, 
	main = 'Classical', xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)')
axis(1, at = p1, labels = 1:10)
# Calculating x*'s x-coordinate
#yardstick = p1[6] - p1[5]
#abline(v = p1[5] + 0.6 * yardstick, lty = 3)

p2 = barplot(100*bpi, 
	main = 'Biased-Coin', xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)')
axis(1, at = p2, labels = 1:10)
# Calculating x*'s x-coordinate
#yardstick = p1[4] - p1[3]
#xfrac = qlogis(btarg, location = exampleMu, scale = exampleSig) - 3
#abline(v = p1[3] + xfrac * yardstick, lty = 3)

dev.off()
