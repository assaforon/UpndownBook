library(cir)
library(upndown)

# Dose-response function used in many examples
exampleMu = 5.6
exampleSig = 2
M = 10

exampleF0 = plogis(1:M, location = exampleMu, scale = exampleSig)
exampleFshift2 = plogis((1:M)-2, location = exampleMu, scale = exampleSig)
exampleFshift4 = plogis((1:M)-4, location = exampleMu, scale = exampleSig)

piF = pivec(exampleF0, matfun = classicmat)
piFs2 = pivec(exampleFshift2, matfun = classicmat)
piFs4 = pivec(exampleFshift4, matfun = classicmat)

pdf('../../output/ch7_boundary_bias.pdf', width = 5, height = 7)
par(mfrow = c(3,1), mar = c(1,1,1,1), mgp = c(2.5, 0.6, 0), tck = -0.01, las = 1, 
				xaxt = 'n', yaxt = 'n', cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)
				
plot(piF, type = 'l', ylab = '', xlab = '', lty = 3, xlim = c(-4, 10) )
points(piF, pch = 19)
abline(v = exampleMu, lty = 2)
abline(v = weighted.mean(1:M, w = piF) )

plot((1:M)-2, piFs2, type = 'l', ylab = '', xlab = '', lty = 3, xlim = c(-4, 10))
points((1:M)-2, piFs2, pch = 19)
abline(v = exampleMu, lty = 2)
abline(v = weighted.mean((1:M)-2, w = piFs2) )

plot((1:M)-4, piFs4, type = 'l', ylab = '', xlab = '', lty = 3, xlim = c(-4, 10))
points((1:M)-4, piFs4, pch = 19)
abline(v = exampleMu, lty = 2)
abline(v = weighted.mean((1:M)-4, w = piFs4) )

dev.off()

