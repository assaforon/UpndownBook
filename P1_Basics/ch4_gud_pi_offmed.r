source('basics_header.r')

wid = 2
g201pi = 100*pivec(exampleF, gudmat, cohort=2, lower=0, upper=1)
g302pi = 100*pivec(exampleF, gudmat, cohort=3, lower=0, upper=2)
g402pi = 100*pivec(exampleF, gudmat, cohort=4, lower=0, upper=2)
g545pi = 100*pivec(exampleF2, gudmat, cohort=5, lower=4, upper=5)
g656pi = 100*pivec(exampleF2, gudmat, cohort=6, lower=5, upper=6)
g767pi = 100*pivec(exampleF2, gudmat, cohort=7, lower=6, upper=7)



#pdf(file.path(outdir, 'ch4_gud_pi_offmed.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

# 30th percentile
plot(g201pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', ylim = c(1, max(g302pi)),
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), lty = 3,
		main = 'Targeting the 30th Percentile')
lines(g302pi, lty = 2, lwd=wid)
lines(g402pi, lty = 1, lwd=wid)

abline(v = qlogis(0.3, location = exampleMu, scale = exampleSig)) 
axis(1, 1:M)


plot(g545pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', ylim = c(1, max(g302pi)), 
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ),
		, main = 'Targeting the 90th Percentile',  lty = 3, xlim = c(5, 16) )
lines(g656pi, lty = 2, lwd=wid)
lines(g767pi, lty = 1, lwd=wid)
#legend('topright', legend = c('(5,2,3)', '(5,1,4)', '(5,0,5)'),lty = c(lines2,'solid'), lwd = 2,
#     bty = 'n', cex = 1.3)
abline(v = qlogis(0.9, location = exampleMu, scale = exampleSig), lty=2, lwd=wid) 

axis(1, 1:(2*M) )

#dev.off()

