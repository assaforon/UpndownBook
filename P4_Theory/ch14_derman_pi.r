source('../P1_Basics/basics_header.r')

btarg1 = 0.3
btarg2 = 0.9

bpi1 = pivec(exampleF, bcdmat, target = btarg1)
bpi2 = pivec(exampleF2, bcdmat, target = btarg2)

derlambda1 = (1 - exampleF[-M])/(exampleF[-1] + 1 - 2*btarg1)
derlambda2 = (2*btarg2 - exampleF2[-(2*M)])/exampleF2[-1]

derpi1 = cumprod(c(1, derlambda1) )
derpi1 = derpi1/sum(derpi1)
derpi2 = cumprod(c(1, derlambda2) )
derpi2 = derpi2/sum(derpi2)

# stop('check!')

pdf(file.path(outdir, 'ch14_derman_pi.pdf'), width = 14, height = 7.3)

par(mfrow = 1:2, mar = c(4,4,4,2), mgp = c(2.5, 0.6, 0), tck = -0.01, 
				las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)
				
lwid = 1.5

plot(100*bpi1, type = 'l', lwd = lwid, xaxt = 'n',
	main = expression(paste(Gamma, ' = 0.3') ), xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)')
axis(1, at =  1:M)	
lines(100*derpi1, lwd = lwid, lty = 2)
abline(v = qlogis(btarg1, location = exampleMu, scale = exampleSig), lty = 3)
legend('topright', lty = 1:2, legend = c('Standard', 'Derman'), bty = 'n', cex = 1.2, lwd = lwid) 

plot(100*bpi2, type = 'l', lwd = lwid, xlim = c(6,15), xaxt = 'n',
	main = expression(paste(Gamma, ' = 0.9') ), xlab = 'Dose-Level', ylab = 'Asymptotic Probability/Frequency (%)')
axis(1, at =  1:(2*M) )	
lines(100*derpi2, lwd = lwid, lty = 2)
abline(v = qlogis(btarg2, location = exampleMu, scale = exampleSig), lty = 3)

dev.off()




