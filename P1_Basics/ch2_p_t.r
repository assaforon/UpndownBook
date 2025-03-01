source('basics_header.r')


btarg = 0.3

bpi = pivec(exampleF, bcdmat, target = btarg)

bp10_1 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 1, n=10)
bp20_1 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 1, n=20)
bp30_1 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 1, n=30)

bp10_4 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 4, n=10)
bp20_4 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 4, n=20)
bp30_4 = currentvec(exampleF, matfun = bcdmat, target=0.3, startdose = 4, n=30)

pdf(file.path(outdir, 'ch2_p_t.pdf'), width = 14, height = 7.3)

layout(t(1:2))

par(mar = c(4,4,4,2), mgp = c(2.5, 0.6, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)

plot(100*bp10_1, type = 'l', main = expression(paste('Biased-Coin, starting at ', d[1])),
	xlab = 'Dose-Level', xaxt = 'n', ylab = 'Allocation Probability (%)', ylim = c(1, 29) )
axis(1, at =  1:10)
lines(100*bp20_1, lty = 3)
lines(100*bp30_1, lty = 2)
lines(100*bpi, lwd = 2)
legend('topright', lty = c(1,3,2,1), lwd = c(1,1,1,2), legend = c('p(10)', 'p(20)', 'p(30)', expression(pi)), bty = 'n', cex = 1.5) 

plot(100*bp10_4, type = 'l', main = expression(paste('Biased-Coin, starting at ', d[4])),
	xlab = 'Dose-Level', xaxt = 'n', ylab = 'Allocation Probability (%)', ylim = c(1, 29))
axis(1, at =  1:10)
lines(100*bp20_4, lty = 3)
lines(100*bp30_4, lty = 2)
lines(100*bpi, lwd = 2)

dev.off()
