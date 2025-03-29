source('basics_header.r')

finex = seq(1, 10, 0.1)
fineF = plogis(finex, location = exampleMu, scale = exampleSig)

btarg = 0.3

pdf(file.path(outdir, 'ch2_balancept.pdf'), width = 14, height = 7.3)

layout(t(1:2))
par(stdpar)

plot(finex, fineF, lty = 2, type = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', 
	main = 'Classical', xlab = 'Dose-Level', ylab = 'Transition Probability')
axis(1, 1:10)
axis(2, (1:10)/10 )
lines(finex, 1-fineF, type = 'l', lwd = 2)
abline(h = 0.5, lty = 3)
abline(v = exampleMu, lty = 3)

plot(finex, fineF, lty = 2, type = 'l', lwd = 2, xaxt = 'n', yaxt = 'n', 
	main = 'Biased-Coin', xlab = 'Dose-Level', ylab = 'Transition Probability')
axis(1, 1:10)
axis(2, (1:10)/10 )
lines(finex, (1-fineF) * btarg / (1-btarg), type = 'l', lwd = 2)
abline(h = btarg, lty = 3)
abline(v = qlogis(btarg, location = exampleMu, scale = exampleSig), lty = 3)

dev.off()
