rm(list=ls())
library(cir)
library(upndown)
library(ggplot2)
theme_set(theme_bw(18))
library(data.table)

# For use on your own machine, change 'outdir' to where you want the figures saved
outdir = ('../../output')

#------------ Prep

load(file.path(outdir, 'biasfigs_makedata.RData'))	

M = dim(l50midmid$scenarios)[1]
n1 = dim(l50midmid$responses)[1]
n2 = dim(l50midmid_long$responses)[1]

### Utility for one-stop-shop bias calculation

fcalc <- function(simdat, directbias = TRUE)
{
	m = dim(simdat$scenarios)[1]
	n = dim(simdat$responses)[1]
#	dcuts = cut(simdat$doses[1:n, ], (0:m)-0.5, labels = 1:m)
	tmp = mapply(function(y,x) sapply(split(y,x), mean), x=split(simdat$doses[1:n,], col(simdat$doses[1:n,])), 
							y=split(simdat$responses, col(simdat$responses)) )
							
	dout = sapply(tmp, function(x,mm) {dout = rep(NA, mm); dout[as.integer(names(x))] <- x; dout}, mm=m)
	if(directbias) return(rowMeans(dout - simdat$scenarios, na.rm=TRUE) )
	dout
}

#------------ Process

# setting up the main plotting dataset
crml50bias = data.table(Design = 'CRM', N = n1, Dose = 1:M, Bias = fcalc(crml50midmid) )			
boinl50bias = data.table(Design = 'Interval (BOIN)', N = n1, Dose = 1:M, Bias = fcalc(boinl50midmid) )				
udl50bias = data.table(Design = 'Classical UDD', N = n1, Dose = 1:M, Bias = fcalc(l50midmid) )
biasdat = rbind(udl50bias, crml50bias, boinl50bias)

crml50biasl = data.table(Design = 'CRM', Dose = 1:M, N = n2, Bias = fcalc(crml50midmid_long) )			
boinl50biasl = data.table(Design = 'Interval (BOIN)', N = n2, Dose = 1:M, Bias = fcalc(boinl50midmid_long) )				
udl50biasl = data.table(Design = 'Classical UDD', N = n2, Dose = 1:M, Bias = fcalc(l50midmid_long) )
biasdat = rbind(biasdat, udl50biasl, crml50biasl, boinl50biasl)

p1 <- ggplot(biasdat, aes(Dose, Bias, lty = Design) ) + 
		geom_rect(ymin = -Inf, ymax = Inf, xmin = M/2, xmax = M/2 + 1, fill = 'gray90') +
		geom_line(lwd = 1) + facet_wrap(~factor(paste('n =', N), levels = paste('n =', c(n1,n2)) ) ) +
		geom_hline(yintercept = 0, lwd = 1.5) + scale_x_continuous(breaks = 1:12) + scale_linetype_manual(values = 1:3) + 
		theme(panel.grid.minor = element_blank() ) + labs(x = 'Dose Level', y = 'Empirical Bias of Observed Response Rates') 
		
ggsave(p1, file = file.path(outdir, 'ch6_biasfig_1.pdf'), width = 12, height = 6)
