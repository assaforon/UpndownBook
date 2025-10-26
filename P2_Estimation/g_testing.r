rm(list=ls())
library(cir)
library(upndown)
outdir = '../../output'
# Overlap CI standard
ovconf = 0.84



#-------------------------------------- Benhamou et al. 2003

# For brevity, we initially use integers to denote the doses. 

xropi03 = c(11:9,10:8,9,10,9,10:7,8:11,10:12,11:7,8,7:10,9,8,9,8:10,9,10,9,10)
xlevo03 = c(11,10,11,10,11:9,10:7,8,7,8:5,6:8,7,8:6,7,6,7,6,7:5,6,7,6:12)

yropi03=(1-diff(xropi03))/2
ylevo03 = (1-diff(xlevo03))/2

xropi03 = xropi03/100
xlevo03 = xlevo03/100

## 84% CI overlap method
ropi03est84 = udest(x = xropi03, y = yropi03, target = 0.5, conf = ovconf, allow1extra = TRUE)
levo03est84 = udest(x = xlevo03, y = ylevo03, target = 0.5, conf = ovconf, allow1extra = TRUE)

## 95% CI of difference method

# ropi03est95 = udest(x = xropi03, y = yropi03, target = 0.5, conf = 0.95, allow1extra = TRUE)
# levo03est95 = udest(x = xlevo03, y = ylevo03, target = 0.5, conf = 0.95, allow1extra = TRUE)

# ropi03widths = c(ropi03est95$point - ropi03est95$lower95conf, ropi03est95$upper95conf - ropi03est95$point)
# levo03widths = c(levo03est95$point - levo03est95$lower95conf, levo03est95$upper95conf - levo03est95$point)
# pointdiff = ropi03est95$point - levo03est95$point

# diff95ci = data.frame(point = pointdiff, lower95conf = pointdiff - sqrt(ropi03widths[1]^2 + levo03widths[2]^2), 
							# upper95conf = pointdiff + sqrt(ropi03widths[2]^2 + levo03widths[1]^2) )
	
diff03ci95 = armDiffCI(xropi03[-40], yropi03, xlevo03[-40], ylevo03, target = 0.5)	

### Plotting

pdf(file.path(outdir, 'Benhamou03.pdf'), width = 12, height = 8)
layout(matrix(1:4, nrow = 2), widths = 3:2 )
par(mar = c(4,4,4,2), mgp = c(2.6, 0.5, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)

dosename1 = "Concentration (%)"
dosevals1 = (5:12) / 100
arms1 = c("Ropivacaine", "Levobupivacaine")
types = c('Trajectory', 'Dose-Response')
udplot(x = xropi03[1:39], y = yropi03, main = paste(arms1[1],types[1]), ytitle = dosename1, doselabels = dosevals1, 
		ylim = range(dosevals1), cex = 1.5 )
legend('bottomright',legend=c('Effective','Ineffective'),pch=c(19,1),bty='n', cex = 1.4)
udplot(x = xlevo03[1:39], y = ylevo03, main = paste(arms1[2],types[1]), ytitle = dosename1, doselabels = dosevals1, 
	ylim = range(dosevals1), cex = 1.5)

drplot(x = xropi03[1:39], y = yropi03, main = paste(arms1[1],types[2]), xtitle = dosename1, ytitle = 'Proportion Effective', doselabels = dosevals1, obsize = .8,
			addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5, curvecol = 1, estcol = 1, curvetype = 2)
drplot(x = xlevo03[1:39], y = ylevo03, main = paste(arms1[2],types[2]), xtitle = dosename1, ytitle = 'Proportion Effective', doselabels = dosevals1,  obsize = .8,
			addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5, curvecol = 1, estcol = 1, curvetype = 2)
dev.off()			
			
#-------------------------------------- Camorcia et al. 2011  (revisited from Averaging chpater)	
	
c11m_x = c(4:6, 5, 4:7, 6, 5, 6:4, 5:8, 7, 8 , 7, 8, 7:9, 8:11, 10, 9)
c11m_y = c( (1 - sign(diff(c11m_x)) )/ 2, 0) 

c11f_x = c(4, 5, 4:6, 5:3, 4:7, rep(6:5, 3), 6, 7:4, 5, 6:4, 5, 6, 5)
c11f_y = c( (1 - sign(diff(c11f_x)) )/ 2, 1) 

c11c_x = c(4:2, rep(3:4, 3), 3, 2:4, 3, 2, 3, 2:3, rep(c(4, 3:5), 3) )
c11c_y = c( (1 - sign(diff(c11c_x)) )/ 2, 1) 

### Testing with "manual FDR":
alpha = 0.05/3

diff11ci95_fc = armDiffCI(c11f_x, c11f_y, c11c_x, c11c_y, target = 0.5, conf = 1 - alpha)	
diff11ci95_mc = armDiffCI(c11m_x, c11m_y, c11c_x, c11c_y, target = 0.5, conf = 1 - 2*alpha)	
diff11ci95_mf = armDiffCI(c11m_x, c11m_y, c11f_x, c11f_y, target = 0.5, conf = 1 - 3*alpha)	

### Plotting

pdf(file.path(outdir, 'Camorcia11_2.pdf'), width = 11, height = 11)
layout(matrix(1:6, nrow = 3), widths = 3:2 )
par(mar = c(4,4,4,2), mgp = c(2.6, 0.5, 0), tck = -0.01, las = 1, cex.lab = 1.8, cex.axis = 1.4, cex.main = 2)

dosename2 = "Bupivacaine (mg)"
dosevals2 = 2:11
arms2 = c("Women (non-pregnant)", "Women (C-Section)", "Men")

udplot(x = c11f_x, y = c11f_y, main = paste(arms2[1],types[1]), ytitle = dosename2, doselabels = dosevals2, 
	ylim = range(dosevals2), cex = 2 )
udplot(x = c11c_x, y = c11c_y, main = paste(arms2[2],types[1]), ytitle = dosename2, doselabels = dosevals2, 
	ylim = range(dosevals2), cex = 2 )
legend('topright',legend=c('Effective','Ineffective'),pch=c(19,1),bty='n', cex = 1.8)
udplot(x = c11m_x, y = c11m_y, main = paste(arms2[3],types[1]), ytitle = dosename2, doselabels = dosevals2, 
	ylim = range(dosevals2), cex = 2 )

drplot(x = c11f_x, y = c11f_y, main = paste(arms2[1],types[2]), xtitle = dosename2, ytitle = 'Proportion Effective',  obsize = 1.,
		doselabels = dosevals2, addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5, curvecol = 1, estcol = 1, curvetype = 2)
drplot(x = c11c_x, y = c11c_y, main = paste(arms2[2],types[2]), xtitle = dosename2, ytitle = 'Proportion Effective', obsize = 1.,
		doselabels = dosevals2, addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5, curvecol = 1, estcol = 1, curvetype = 2)
drplot(x = c11m_x, y = c11m_y, main = paste(arms2[3],types[2]), xtitle = dosename2, ytitle = 'Proportion Effective', obsize = 1.,
		doselabels = dosevals2, addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5, curvecol = 1, estcol = 1, curvetype = 2)
dev.off()		
