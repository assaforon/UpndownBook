rm(list=ls())
library(cir)
library(upndown)

# For use on your own machine, change 'outdir' to where you want the figures saved
outdir = ('../../output')

#  **An up-and-down experiment that has generated some controversy**
#  
# Van Elstraete, AC et al. The Median Effective Dose of Preemptive Gabapentin 
#      on Postoperative Morphine Consumption After Posterior Lumbar Spinal Fusion. 
#      *Anesthesia & Analgesia* 2008, 106: 305-308.

# It was a classical median-finding up-and-down study.

van08x = c(4:7, 6:13, 12:19, 18:21, 20, 19:23, 22, 21:23, 22:19, 20:23,
          22:24, 23, 22, 23, 22:25, 24:22, rep(23:24,2), 23, 22)
# With U&D, responses (except the last one) can be read off the doses:
van08y = c( (1 - sign(diff(van08x)))/2, 0 )

# A doseResponse object binding the x and y together is helpful
#     with many cir functions

van08dr = doseResponse(x=van08x, y=van08y)



pdf(file.path(outdir,'isotonic1.pdf'),width = 6, height = 8)

par(mfrow=2:1, mar=c(4,4,1,1), mgp=c(2.5,0.8,0), cex.axis = 0.8, las = 1)

udplot(van08x, van08y, xtitle = "Patient Number", ytitle = 'Gabapentin (mg/kg)' )

drplot(van08x, van08y, addcurve = TRUE, addest = TRUE, target = 0.5,
		percents = TRUE, xtitle = 'Gabapentin (mg/kg)', ytitle = "Percent Effective")

van08IR = oldPAVA(van08dr, full = TRUE)
lines(I(100*y) ~ x, data = van08IR$output, lty = 2)
	
dev.off()

###### Now using this study to visualize interval estimates

van08CIR = cirPAVA(van08dr,full = TRUE, adaptiveShrink = TRUE, target = 0.5)

van08fwd = quickIsotone(van08dr, adaptiveShrink = TRUE, target = 0.5, 
				outx = van08CIR$shrinkage$x)
van08inv = quickInverse(van08dr, adaptiveShrink = TRUE, target = 0.5) 

pdf(file.path(outdir,'isotonic5.pdf'),width = 6, height = 5)

par(mar=c(4,4,1,1), mgp=c(2.5,0.8,0), cex.axis = 0.8, las = 1)

plot(y ~ x, data = van08fwd, type = 'l', ylim = 0:1,
      xlab = 'Gabapentin (mg/kg)', ylab = "Proportion Effective")
lines(lower90conf ~ x, data = van08fwd, lty = 2)
lines(upper90conf ~ x, data = van08fwd, lty = 2)
points(target ~ point, data = van08inv, pch = 19, cex = 2)
lines(c(van08inv$lower90conf, van08inv$upper90conf), c(.5, .5), lwd = 1.5)

dev.off()

