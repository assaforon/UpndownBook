source('basics_header.r')

wid = 2
g201pi = 100*pivec(exampleF, gudmat, cohort=2, lower=0, upper=1)
g323pi = 100*pivec(exampleF2, gudmat, cohort=3, lower=2, upper=3)

k12pi = 100*pivec(exampleF, kmatMarg, k=2, lowTarget = TRUE)
k31pi = 100*pivec(exampleF2, kmatMarg, k=3, lowTarget = FALSE)


pdf(file.path(outdir, 'ch4_kr_gud_pi.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

# 30th percentile
plot(g201pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', ylim = c(1, max(g201pi)),
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), lty = 2,
		main = expression(paste(KRD[1/2],'     (F* = 0.293)') ), xlim = c(1, 9) )
lines(k12pi, lty = 1, lwd=wid)
legend('topright', legend = c('KRD', 'GUD / Arrivals', 'x*'),lty = 1:3, lwd = c(wid, wid, 1),
     bty = 'n', cex = 1.3)

abline(v = qlogis(k2targ(2, lowTarget = TRUE), location = exampleMu, scale = exampleSig), lty=3) 
axis(1, 1:M)


plot(g323pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', ylim = c(1, max(g323pi)), 
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ),
		, main = expression(paste(KRD[3/1],'     (F* = 0.794)') ),  lty = 2, xlim = c(3, 14) )
lines(k31pi, lty = 1, lwd=wid)
legend('topright', legend = c('KRD', 'GUD / Arrivals', 'x*'),lty = 1:3, lwd = c(wid, wid, 1),
     bty = 'n', cex = 1.3)
abline(v = qlogis(k2targ(3, lowTarget = FALSE), location = exampleMu, scale = exampleSig), lty=3) 

axis(1, 1:(2*M) )

dev.off()

