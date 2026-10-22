cat(base::date(), '\n')
rm(list=ls())
library(upndown)
outdir = '../output'

chosen = 15

wid = 12
hgt = 10

load(file.path(outdir, 'grandsim90w.RData'))	


layout(t(1:2)) #, widths = 3:2)
forecol = 'white'
scol = 'yellow'

par(las = 1, tck = -0.01, bg = '#004033', fg = forecol, col.axis = forecol, col.lab = forecol,
	col.main = forecol, cex.lab = 1.6, cex.axis = 1.2) #, mgp = c(2.5,0.7,0) )

udplot(kw90himid$doses[-61,chosen], kw90himid$responses[ ,chosen], cex = 1.5, xaxt = 'n', symbcol = forecol)
axis(1, at = seq(5, 60, 5) )
drplot(kw90himid$doses[-61,chosen], kw90himid$responses[ ,chosen], target = .9, balancePt = k2targ(6, lowTarget = FALSE), 
			addcurve = TRUE, addest = TRUE, curvethick = 3, estcol = 'hotpink', obsize = 1, esthick = 2,
			ytitle = 'Positive Response Rate (%)', curvecol = forecol, symbcol = scol, percents = TRUE)


