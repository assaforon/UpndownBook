rm(list=ls())
library(cir)
library(upndown)
outdir = '../../output'


# For brevity, we initially use integers to denote the doses. 

xropi03 = c(11:9,10:8,9,10,9,10:7,8:11,10:12,11:7,8,7:10,9,8,9,8:10,9,10,9,10)
xlevo03 = c(11,10,11,10,11:9,10:7,8,7,8:5,6:8,7,8:6,7,6,7,6,7:5,6,7,6:12)

yropi03=(1-diff(xropi03))/2
ylevo03 = (1-diff(xlevo03))/2

xropi03 = xropi03/100
xlevo03 = xlevo03/100

## 84% CI overlap method
ropi03est84 = udest(x = xropi03, y = yropi03, target = 0.5, conf = 0.84, allow1extra = TRUE)
levo03est84 = udest(x = xlevo03, y = ylevo03, target = 0.5, conf = 0.84, allow1extra = TRUE)

## 95% CI of difference method

ropi03est95 = udest(x = xropi03, y = yropi03, target = 0.5, conf = 0.95, allow1extra = TRUE)
levo03est95 = udest(x = xlevo03, y = ylevo03, target = 0.5, conf = 0.95, allow1extra = TRUE)

ropi03widths = c(ropi03est95$point - ropi03est95$lower95conf, ropi03est95$upper95conf - ropi03est95$point)
levo03widths = c(levo03est95$point - levo03est95$lower95conf, levo03est95$upper95conf - levo03est95$point)
pointdiff = ropi03est95$point - levo03est95$point

diff95ci = data.frame(point = pointdiff, lower95conf = pointdiff - sqrt(ropi03widths[1]^2 + levo03widths[2]^2), 
							upper95conf = pointdiff + sqrt(ropi03widths[2]^2 + levo03widths[1]^2) )
							
### Plotting

par(mfrow = 1:2, mar = c(4,4,4,2), mgp = c(2.6, 0.5, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)
dosename = "Concentration (%)"
dosevals = (5:12) / 100

udplot(x = xropi03[1:39], y = yropi03, main = "Ropivacaine", ytitle = dosename, doselabels = dosevals, ylim = range(dosevals) )
legend('bottomright',legend=c('Effective','Ineffective'),pch=c(19,1),bty='n', cex = 1.2)
udplot(x = xlevo03[1:39], y = ylevo03, main = "Levobupivacaine", ytitle = dosename, doselabels = dosevals, ylim = range(dosevals))