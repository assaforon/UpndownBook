source('basics_header.r')

pdf(file.path(outdir, 'ch1_doseresponse.pdf'), width = 10, height = 7)
par(mar = c(4,4,1,1), mgp = c(2.5, 0.6, 0), tck = -0.01)

expandedX = seq(-3, 15, 0.1)

plot(expandedX, plogis(expandedX, location = exampleMu, scale = exampleSig), type = 'n', xaxt = 'n', 
				las = 1, xlab = 'x', ylab = 'Pr (Y = 1)', cex.lab = 1.7, cex.axis = 1.2, bty = 'l' )
				
# Parameters for circles

csize = 12
cchar = 21
				
points(qlogis(0.3, location = exampleMu, scale = exampleSig), 0.3, pch = cchar, cex = csize, col = 'white', bg = 'gray90')
points(qlogis(0.8, location = exampleMu, scale = exampleSig), 0.8, pch = cchar, cex = csize, col = 'white', bg = 'gray40')

points(exampleMu, 0.5, pch = cchar, cex = csize, col = 'white', bg = 'gray70')



lines(expandedX, plogis(expandedX, location = exampleMu, scale = exampleSig), lwd = 3)

text(-0., 0.09, label = 'F(x)', cex = 1.5, srt = 15)

dev.off()