source('simother_header.r')
library(upndown)
library(ggplot2)
library(ggtext)
theme_set(theme_bw(16)) 
#

outdir = 'C:/GitHub/output'
mycurve = 1
M = 8

weib30parm = fread(file.path(outdir, 'scenarios_weib30.csv') )

weib30F0 = weib30parm[1:10,][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)

cum30 = cumulvec(weib30F0[,mycurve], matfun = kmatMarg, k = 2, lowTarget = TRUE, startdose = 1, n = 30)
vec30 = currentvec(weib30F0[,mycurve], matfun = kmatMarg, k = 2, lowTarget = TRUE, startdose = 1, n = 30)
vec15 = currentvec(weib30F0[,mycurve], matfun = kmatMarg, k = 2, lowTarget = TRUE, startdose = 1, n = 16)
vecpi = pivec(weib30F0[,mycurve], matfun = kmatMarg, k = 2, lowTarget = TRUE)

dists = data.frame(dose=1:M, n30 = 100*cum30, v15 = 100*vec15, v30 = 100*vec30, 
		asym = 100*vecpi, dr = 100*weib30F0[,mycurve])

ggplot(dists, aes(x=dose)) + geom_line(aes(y=dr), lwd = 5, col = 'green', alpha = 0.25) + 
		scale_x_continuous(breaks = 1:8) + geom_col(aes(y=n30),fill='steelblue2', width=0.5) -> p1
p2 <- p1 +  geom_point(aes(y=v15),col='red', size=4) + geom_line(aes(y=v15),col='red', lty=2) + 
	geom_point(aes(x=dose, y=asym), size=4) +  geom_line(aes(y=asym), lty=2) + 
	coord_cartesian(ylim = c(0,43), expand = 0) + 
	labs(x = "Dose-Level (arbitrary units)", y = "Allocation Proportion (%)") +
# , "<span style = 'font-family: Arial; font-size: 12pt'>[or response rate]</span>") ) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(p2, file = file.path(outdir, 'NEJarticle_alloc.pdf'), width = 8, height = 8) 
# stop('vec')
#--------------------- CIR plot

dlabel = expression(paste('Phenylephrine dose (', mu, 'g)' ) )
george10x = 80 + 20 * c(1, rep(2, 5), 1, 1, 0, 0, rep(1, 7), 0:2, 2, 2, rep(1, 4), 2, 1, 1, 2, 2, 
                        rep(3, 5), 4, 5, 5, rep(4, 6))
george10y = c(ifelse(diff(george10x) > 0, 0, 1), 1)

pdf(file.path(outdir, 'NEJarticle_CIR.pdf'), width = 8, height = 7.5)

par(las = 1, cex.axis = 1.3, cex.lab = 1.6, mar = c(4.5,4.5,1,1) )
drplot(x=george10x, y=george10y, addest = TRUE, target = 0.9, addcurve = TRUE, 
			balancePt = 10/11, xtitle = dlabel, percents = TRUE, ytitle = "Efficacy (%)" )

irest = quickIsotone(x=george10x, y=george10y, estfun = oldPAVA)

lines((100*y)~x, data=irest, lty=2)

legend('bottomright', bty = 'n', lty = c(0,2,1,1), pch =  c(4,NA,NA,19), col = c(1,1,'blue','purple'),
				legend = c("Observed Rates", 'IR Curve', 'CIR Curve', 'Target Estimate'), cex = 1.4, lwd = 2)
				
dev.off()

