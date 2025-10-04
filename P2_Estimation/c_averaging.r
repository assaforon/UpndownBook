rm(list=ls())
library(upndown)

# For use on your own machine, change 'outdir' to where you want the figures saved
outdir = ('../../output')

 
# =============================================

#### Reversals and estimates

# =============================================

## Plotting constants

grey = 'grey60'
ssize = 1.5

# Camorcia et al. 2011
# "Effect of sex and pregnancy on the potency of intrathecal bupivacaine..."
# European Journal of Anaesthesiology 2011, Volume 28, p. 240-244

# There were two female arms (c-section and non-pregnant) 
#   and one male arm

c11m_x = c(4:6, 5, 4:7, 6, 5, 6:4, 5:8, 7, 8 , 7, 8, 7:9, 8:11, 10, 9)
c11m_xplus = c(c11m_x, 10)
c11m_y = c( (1 - sign(diff(c11m_x)) )/ 2, 0) 
c11m_rev = reversmean(c11m_xplus, c11m_y, conf = NULL)
c11m_weth = reversmean(c11m_xplus, c11m_y, all = FALSE, rstart = 1, conf = NULL)
c11m_dyn = dynamean(c11m_xplus, conf = NULL)

c11f_x = c(4, 5, 4:6, 5:3, 4:7, rep(6:5, 3), 6, 7:4, 5, 6:4, 5, 6, 5)
c11f_y = c( (1 - sign(diff(c11f_x)) )/ 2, 1) 
c11f_xplus = c(c11f_x, 4)
c11f_rev = reversmean(c11f_xplus, c11f_y, conf = NULL)
c11f_weth = reversmean(c11f_xplus, c11f_y, all = FALSE, rstart = 1, conf = NULL)
c11f_dyn = dynamean(c11f_xplus, conf = NULL)

c11c_x = c(4:2, rep(3:4, 3), 3, 2:4, 3, 2, 3, 2:3, rep(c(4, 3:5), 3) )
c11c_xplus = c(c11c_x, 4)
c11c_y = c( (1 - sign(diff(c11c_x)) )/ 2, 1) 
c11c_rev = reversmean(c11c_xplus, c11c_y, conf = NULL)
c11c_weth = reversmean(c11c_xplus, c11c_y, all = FALSE, rstart = 1, conf = NULL)
c11c_dyn = dynamean(c11c_xplus, conf = NULL)

#------------- Plots

pdf(file.path(outdir,'averaging2.pdf'),width = 8, height = 7)
layout(1:3, heights = c(6,6,7)) # Lowest frame bigger for x-title

par(tcl=-0.4, mar=c(1,4,2,1), mgp=c(2.5,1.,0), las=1, cex.lab = 1.4, cex.main = 1.5)

# 1
udplot(c11f_x, c11f_y, shape='square', connect = FALSE, 
	main = "Women (non-pregnant)", xtitle = "", ytitle = 'Bupivacaine (mg)', 
	cex = ssize, ylim = c(2,11), xaxt = 'n', doselabels = 2:11) # xlim = c(1, 31), 
points(31, c11f_xplus[31], pch = 22, col = grey, bg = grey, cex = ssize)
# Showing reversals 
frevs = reversals(x = c11f_x, y = c11f_y)
points(frevs, c11f_x[frevs], cex = 3.5)

# Estimates
abline(h = c11f_rev)
abline(h = c11f_weth, lty = 3)
abline(h = c11f_dyn, lty = 2, lwd = 1.5)

# 2
udplot(c11c_x, c11c_y, shape='square', connect = FALSE,
 main = "Women (C-Section)", xtitle = "", ytitle = 'Bupivacaine (mg)',
 cex = ssize, ylim = c(2,11), xaxt = 'n', doselabels = 2:11) # xlim = c(1, 31), 
points(31, c11c_xplus[31], pch = 22, col = grey, bg = grey, cex = ssize)
crevs = reversals(x = c11c_x, y = c11c_y)
points(crevs, c11c_x[crevs], cex = 3.5)

# Estimates
abline(h = c11c_rev)
abline(h = c11c_weth, lty = 3)
abline(h = c11c_dyn, lty = 2, lwd = 1.5)

# 3
par(mar=c(4,4,2,1))
udplot(c11m_x, c11m_y, shape='square', connect = FALSE,
 main = "Men", xtitle = "Patient Number", ytitle = 'Bupivacaine (mg)',
 cex = ssize, ylim = c(2,11), doselabels = 2:11) # xlim = c(1, 31), 
points(31, c11m_xplus[31], pch = 22, col = grey, bg = grey, cex = ssize)
mrevs = reversals(x = c11m_x, y = c11m_y)
points(mrevs, c11m_x[mrevs], cex = 3.5)

# Estimates
abline(h = c11m_rev)
abline(h = c11m_weth, lty = 3)
abline(h = c11m_dyn, lty = 2, lwd = 1.5)

dev.off()


### Demonstrating the dynamic estimator

# Rerunning the estimates w/full = TRUE

c11f_dyn = dynamean(c11f_xplus, full = TRUE, conf = NULL)
c11m_dyn = dynamean(c11m_xplus, full = TRUE, conf = NULL)

pdf(file.path(outdir,'averaging3.pdf'),width = 10, height = 8)
layout(1:2, heights = c(6,7)) # Lowest frame bigger for x-title

par(tcl=-0.4, mar=c(1,4,2,1), mgp=c(2.5,1.,0), las=1, cex.lab = 1.4, cex.main = 1.5)


udplot(c11f_x, c11f_y, shape='square', connect = FALSE, 
 xtitle = "", ytitle = 'Bupivacaine (mg)', main = "Women (non-pregnant)",
 cex = ssize, ylim = c(2,11), xaxt = 'n', doselabels = 2:11) # xlim = c(1, 31), 
points(31, c11f_xplus[31], pch = 22, col = grey, bg = grey, cex = ssize)
for(a in 1:c11f_dyn$startpt) points(a, c11f_dyn$signsmeans[1,a+1], pch = '_', cex=2)
abline(v = 14.5, lwd = 1.5, lty = 3)

par(mar=c(4,4,2,1))

udplot(c11m_x, c11m_y, shape='square', connect = FALSE, 
 xtitle = "Patient Number", ytitle = 'Bupivacaine (mg)',  main = "Men",
 cex = ssize, ylim = c(2,11), doselabels = 2:11) # xlim = c(1, 31), 
points(31, c11m_xplus[31], pch = 22, col = grey, bg = grey, cex = ssize)
for(a in 1:c11m_dyn$startpt) points(a, c11m_dyn$signsmeans[1,a+1], pch = '_', cex=2)
abline(v = 14.5, lwd = 1.5, lty = 3)

dev.off()


#------- Last but not least: Bootstrap F

library(cir)


c11m_dr = doseResponse(x=c11m_x, y=c11m_y)

drplot(c11m_x, c11m_y, target = 0.5,
		xtitle = 'Bupivacaine (mg)', ytitle = "Proportion Effective")
		
c11m_IR = oldPAVA(c11m_dr)

c11m_CIR = cirPAVA(c11m_dr,full = TRUE, adaptiveShrink = TRUE, target = 0.5)

c11m_boot = cirPAVA(c11m_dr,full = TRUE, adaptiveShrink = TRUE, target = 0.5, nmin = 1)

xrange = range(c11m_boot$output$x)
yrange = range(c11m_boot$output$y)

yadd = c(yrange[1]/2, (1+yrange[2])/2)


pdf(file.path(outdir,'averaging4.pdf'),width = 8, height = 4.5)

par(tcl=-0.4, mar=c(4,4,1,1), mgp=c(2.5,1.,0), las=1, cex.lab = 1.2)

drplot(c11m_x, c11m_y,
		xtitle = 'Bupivacaine (mg)', ytitle = "Proportion Effective", xlim = c(xrange[1]-2, xrange[2]+2), las = 1, xaxt = 'n' ) 
axis(1, at = 2:13)

lines(y~x, data = c11m_CIR$shrinkage)		
lines(c11m_dr$x, c11m_IR, lty = 3, lwd = 1.5)

lines(c(xrange[1]-2, xrange[1]-1, c11m_boot$output$x, xrange[2]+1, xrange[2]+2),
     c(0, yadd[1], c11m_boot$output$y, yadd[2], 1), lty = 2, lwd = 2)
	 
legend('topleft', legend = c('IR (naive)', 'CIR (standard)', 'CIR (for bootstrap)'), lty = c(3,1,2), bty = 'n')

dev.off()




		
		





