cat(base::date(), '\n')
rm(list=ls())
source('simulation_header.r')
library(magrittr)

M = 12

grays = c('black', paste('gray', seq(10,95,5), sep='') )
ncolor = length(grays)
ncurves = 50

weib50parm = fread(file.path(outdir, 'scenarios_weib50.csv') )
logi50parm = fread(file.path(outdir, 'scenarios_logi50.csv') )


weib50F = weib50parm[1:ncurves, ][ ,as.vector(pweib3(1:M, shap, scal, shif+offs)), 
					by = 'row0'] %$% matrix(V1, nrow = M)
logi50F = logi50parm[1:ncurves, ][ ,as.vector(plogis(1:M, location=loc, scale=scal)), 
					by = 'row0'] %$% matrix(V1, nrow = M)


pdf(file.path(outdir, 'ensemble_plots.pdf'), width = 9, height = 5)

layout(t(1:2), widths = 20:19)
par(mar = c(3, 3., 3, 1), mgp = c(1.7, 0.4, 0), las = 1, tck = -0.01)

plot(c(1,M),0:1, type='n', xlab = 'Dose Level', ylab = 'F', main = 'Weibull')
for(a in 1:ncurves) lines(1:M, weib50F[ , a], col= grays[a %% ncolor + 1])
abline(h=0.5, lty=2)

par(mar = c(3, 1, 3, 1) )

plot(c(1,M),0:1, type='n', xlab = 'Dose Level', ylab = '', main = 'Logistic')
for(a in 1:ncurves) lines(1:M, logi50F[ , a], col= grays[a %% ncolor + 1])
abline(h=0.5, lty=2)

dev.off()
cat(base::date(), '\n')
