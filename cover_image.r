cat(base::date(), '\n')
rm(list=ls())
library(Hmisc)
library(upndown)
outdir = '../output'

chosen = 15

wid = 12
hgt = 10

load(file.path(outdir, 'grandsim90w.RData'))	

ruler = 4
png('../output/cover_image.png', width = 1100*ruler, height = 500*ruler)

layout(t(1:2), widths = 9:7)
forecol = 'white'
scol = 'yellow'

par(las = 1, tck = -0.005, bg = '#003228', fg = forecol, col.axis = forecol, col.lab = forecol, bty = 'l', # lend = 1,
	col.main = forecol, cex.lab = 1.4*ruler, cex.axis = 1.2*ruler, mar = c(15,16,3,3) , mgp = c(10,2.6,0) )

udplot(kw90himid$doses[-61,chosen], kw90himid$responses[ ,chosen], symbsize = 1.8*ruler, xaxt = 'n', 
			symbcol = scol, linecol = forecol, obslabels = 100, lwd = 1.*ruler,  )

mgp.axis.labels(c(5,2.,0), type='x') 
mgp.axis(1, at = seq(5, 60, 5) )

drplot(kw90himid$doses[-61,chosen], kw90himid$responses[ ,chosen], target = .9, balancePt = k2targ(6, lowTarget = FALSE), 
			addcurve = TRUE, addest = TRUE, curvethick = 2*ruler, estcol = 'magenta', obsize = 1*ruler, 
			esthick = 2*ruler, estsize = 4*ruler, xaxt = 'n',
			ytitle = 'Positive Response Rate (%)', curvecol = forecol, symbcol = scol, percents = TRUE)

mgp.axis(1, at = 4:10 )

dev.off()
