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

nmin = 1 # Minimum n for shrinkage fix
n_s = 1 # shrinkage weight

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

## Similarly, extracting N_m

ncalc <- function(simdat)
{
	m = dim(simdat$scenarios)[1]
	n = dim(simdat$responses)[1]
#	dcuts = cut(simdat$doses[1:n, ], (0:m)-0.5, labels = 1:m)
	tmp = apply(simdat$doses[1:n,], 2, table)							
	sapply(tmp, function(x,mm) {dout = rep(NA, mm); dout[as.integer(names(x))] <- x; dout}, mm=m)
}

#------------ Process

### setting up the main plotting dataset for p1
crml50bias = data.table(Design = 'CRM', N = n1, Dose = 1:M, Bias = fcalc(crml50midmid) )			
boinl50bias = data.table(Design = 'Interval (BOIN)', N = n1, Dose = 1:M, Bias = fcalc(boinl50midmid) )				
udl50bias = data.table(Design = 'Classical UDD', N = n1, Dose = 1:M, Bias = fcalc(l50midmid) )
biasdat = rbind(udl50bias, crml50bias, boinl50bias)

crml50biasl = data.table(Design = 'CRM', Dose = 1:M, N = n2, Bias = fcalc(crml50midmid_long) )			
boinl50biasl = data.table(Design = 'Interval (BOIN)', N = n2, Dose = 1:M, Bias = fcalc(boinl50midmid_long) )				
udl50biasl = data.table(Design = 'Classical UDD', N = n2, Dose = 1:M, Bias = fcalc(l50midmid_long) )
biasdat = rbind(biasdat, udl50biasl, crml50biasl, boinl50biasl)

### p2 is with bias mitigation. Setting it up is a tad more complicated:

# First, extracting R_m and N_M
udl50temp = fcalc(l50midmid, directbias = FALSE)
udl50n = ncalc(l50midmid)
crml50temp = fcalc(crml50midmid, directbias = FALSE)
crml50n = ncalc(crml50midmid)
boinl50temp = fcalc(boinl50midmid, directbias = FALSE)
boinl50n = ncalc(boinl50midmid)

udl50temp_long = fcalc(l50midmid_long, directbias = FALSE)
udl50n_long = ncalc(l50midmid_long)
crml50temp_long = fcalc(crml50midmid_long, directbias = FALSE)
crml50n_long = ncalc(crml50midmid_long)
boinl50temp_long = fcalc(boinl50midmid_long, directbias = FALSE)
boinl50n_long = ncalc(boinl50midmid_long)

# Now, applying the mitigation

udl50fix = ifelse(udl50n > nmin, (udl50temp*udl50n + n_s * 0.5)/(udl50n + n_s), udl50temp)
crml50fix = ifelse(crml50n > nmin, (crml50temp*crml50n + n_s * 0.5)/(crml50n + n_s), crml50temp)
boinl50fix = ifelse(boinl50n > nmin, (boinl50temp*boinl50n + n_s * 0.5)/(boinl50n + n_s), boinl50temp)
udl50fix_long = ifelse(udl50n_long > nmin, (udl50temp_long*udl50n_long + n_s * 0.5)/(udl50n_long + n_s), udl50temp_long)
crml50fix_long = ifelse(crml50n_long > nmin, (crml50temp_long*crml50n_long + n_s * 0.5)/(crml50n_long + n_s), crml50temp_long)
boinl50fix_long = ifelse(boinl50n_long > nmin, (boinl50temp_long*boinl50n_long + n_s * 0.5)/(boinl50n_long + n_s), boinl50temp_long)

stopifnot(all(l50midmid$scenarios == crml50midmid$scenarios))
stopifnot(all(l50midmid$scenarios == l50midmid_long$scenarios))
stopifnot(all(l50midmid_long$scenarios == crml50midmid_long$scenarios))

# Finally, setting up the fixed data for plotting
# rowMeans(dout - l50midmid$scenarios, na.rm=TRUE) )

udl50fias = data.table(Design = 'Classical UDD', N = n1, Dose = 1:M, 
		Bias = rowMeans(udl50fix - l50midmid$scenarios, na.rm=TRUE) )		
boinl50fias = data.table(Design = 'Interval (BOIN)', N = n1, Dose = 1:M, 
		Bias = rowMeans(boinl50fix - l50midmid$scenarios, na.rm=TRUE) )				
crml50fias = data.table(Design = 'CRM', N = n1, Dose = 1:M, 		
		Bias = rowMeans(crml50fix - l50midmid$scenarios, na.rm=TRUE) )		
fiasdat = rbind(udl50fias, crml50fias, boinl50fias)

udll50fias = data.table(Design = 'Classical UDD', N = n2, Dose = 1:M, 
		Bias = rowMeans(udl50fix_long - l50midmid$scenarios, na.rm=TRUE) )		
boinll50fias = data.table(Design = 'Interval (BOIN)', N = n2, Dose = 1:M, 
		Bias = rowMeans(boinl50fix_long - l50midmid$scenarios, na.rm=TRUE) )				
crmll50fias = data.table(Design = 'CRM', N = n2, Dose = 1:M, 		
		Bias = rowMeans(crml50fix_long - l50midmid$scenarios, na.rm=TRUE) )		
fiasdat = rbind(fiasdat, udll50fias, crmll50fias, boinll50fias)


#---------------- PLot
p1 <- ggplot(biasdat, aes(Dose, Bias) ) + # geom_point(size=1.3) +
		geom_rect(ymin = -Inf, ymax = Inf, xmin = M/2, xmax = M/2 +1, fill = 'gray90') +
		geom_line(aes(lty = Design), lwd = 1) + facet_wrap(~factor(paste('n =', N), levels = paste('n =', c(n1,n2)) ) ) +
		geom_hline(yintercept = 0, lwd = 1.5) + scale_x_continuous(breaks = 1:12) + scale_linetype_manual(values = 1:3) + 
		theme(panel.grid.minor = element_blank() ) + labs(x = 'Dose Level', y = 'Observed Bias of Response Rates') 
		
ggsave(p1, file = file.path(outdir, 'ch6_biasfig_1.pdf'), width = 12, height = 6)

p2 <- ggplot(fiasdat, aes(Dose, Bias) ) + # geom_point(size=1.3) +
		geom_rect(ymin = -Inf, ymax = Inf, xmin = M/2, xmax = M/2 +1, fill = 'gray90') +
		geom_line(aes(lty = Design), lwd = 1) + facet_wrap(~factor(paste('n =', N), levels = paste('n =', c(n1,n2)) ) ) +
		geom_hline(yintercept = 0, lwd = 1.5) + scale_x_continuous(breaks = 1:12) + scale_linetype_manual(values = 1:3) + 
		theme(panel.grid.minor = element_blank() ) + labs(x = 'Dose Level', y = 'Observed Bias of Response Rates') + 
		ylim(range(biasdat$Bias)) # Keeping same vertical scale
		
ggsave(p2, file = file.path(outdir, 'ch6_biasfig_2_fix.pdf'), width = 12, height = 6)

