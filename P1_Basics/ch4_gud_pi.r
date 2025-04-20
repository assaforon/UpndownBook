source('basics_header.r')

layout(t(1:2))
par(stdpar)

wid = 2
cpi = 100*pivec(exampleF, classicmat)
g202pi = 100*pivec(exampleF, gudmat, cohort=2, lower=0, upper=2)
g303pi = 100*pivec(exampleF, gudmat, cohort=3, lower=0, upper=3)
g404pi = 100*pivec(exampleF, gudmat, cohort=4, lower=0, upper=4)
g505pi = 100*pivec(exampleF, gudmat, cohort=5, lower=0, upper=5)
g514pi = 100*pivec(exampleF, gudmat, cohort=5, lower=1, upper=4)
g523pi = 100*pivec(exampleF, gudmat, cohort=5, lower=2, upper=3)



pdf(file.path(outdir, 'ch4_gud_pi.pdf'), width = 14, height = 7.3)
layout(t(1:2))
par(stdpar)

lines1 = c('19', '27', '66', '93')
lines2 = c('1955', '4694')
# (K,0,K) family
plot(g505pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', 
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ), main = '(K,0,K) GUDs')
lines(cpi, lty = lines1[1], lwd=wid)
lines(g202pi, lty = lines1[2], lwd=wid)
lines(g303pi, lty = lines1[3], lwd=wid)
lines(g404pi, lty = lines1[4], lwd=wid)
legend('topright', legend = c('Classical (K=1)', paste('K =', 2:5)), lty = c(lines1,'solid'), lwd = 2, 
		bty = 'n', cex = 1.3)
axis(1, 1:M)


plot(g505pi, type = 'l', xaxt = 'n', lwd=wid, xlab = 'Dose-Level', 
		ylab = expression(paste('Stationary/Asymptotic Distribution  ',pi, ' (%)' ) ),
		, main = "Median-Targeting K=5 GUDs")
lines(g514pi, lty = lines2[2], lwd=wid)
lines(g523pi, lty = lines2[1], lwd=wid)
legend('topright', legend = c('(5,2,3)', '(5,1,4)', '(5,0,5)'),lty = c(lines2,'solid'), lwd = 2,
     bty = 'n', cex = 1.3)
axis(1, 1:M)

dev.off()

