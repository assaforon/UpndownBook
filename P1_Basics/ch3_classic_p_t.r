source('basics_header.r')


cpi = pivec(exampleF, classicmat)

cp9_1 = currentvec(exampleF, matfun = classicmat, startdose = 2, n=9)
cp10_1 = currentvec(exampleF, matfun = classicmat, startdose = 2, n=10)
cp29_1 = currentvec(exampleF, matfun = classicmat, startdose = 2, n=29)
cp30_1 = currentvec(exampleF, matfun = classicmat, startdose = 2, n=30)

#pdf(file.path(outdir, 'ch3_classic_p_t.pdf'), width = 14, height = 7.3)

layout(t(1:2))

par(mar = c(4,4,4,2), mgp = c(2.5, 0.6, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)

plot(100*cp9_1, type = 'l', lty = 2, main = expression(paste('Classical UDD, starting at ', d[2])),
	xlab = 'Dose-Level', xaxt = 'n', ylab = 'Allocation Probability (%)', ylim = c(1.8, 51) )
axis(1, at =  1:10)
lines(100*cp10_1, lty = 3)
lines(50*(cp10_1 + cp9_1))
lines(100*cpi, lwd = 2)

legend('topleft', lty = c(2,3,1,1), lwd = c(1,1,1,2), legend = c('p(9)', 'p(10)', '[p(9)+p(10)]/2', expression(pi)), bty = 'n', cex = 1.2) 

plot(100*cp29_1, type = 'l', lty = 2, main = expression(paste('Classical UDD, starting at ', d[2])),
	xlab = 'Dose-Level', xaxt = 'n', ylab = 'Allocation Probability (%)', ylim = c(1.8, 51) )
axis(1, at =  1:10)
lines(100*cp30_1, lty = 3)
lines(50*(cp30_1 + cp29_1))
lines(100*cpi, lwd = 2)

legend('topleft', lty = c(2,3,1,1), lwd = c(1,1,1,2), legend = c('p(29)', 'p(30)', '[p(29)+p(30)]/2', expression(pi)), bty = 'n', cex = 1.2) 

#dev.off()
