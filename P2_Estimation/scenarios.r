source('simulation_header.r')
source("C:\\Documents\\dosefinding\\GenericSimCode\\MichelesPalette.r")

nsim = 1000
midtarg = 6.5

weib90 = data.table(scal = runif(nsim, 2, 10),
					# shap = runif(nsim, .5, 10),
					shap = rgamma(nsim, shape = 2.5, scale = 1.5) + 0.1,
					offs = runif(nsim, -.5, .5) )
					
weib90[ , shif := weibshift(shp=shap, scl=scal, targx = midtarg, targy = 0.9) ]					
weib90[ , t90 := qweib3(0.9, shp=shap, scl=scal, shift = shif+offs) ]					

weib90[ , pp1 := pweib3(midtarg+0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib90[ , pm1 := pweib3(midtarg-0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib90[ , pm2 := pweib3(midtarg-1.5, shp=shap, scl=scal, shift = shif+offs) ]

# plot(c(1,10),0:1, type='n')
# for(a in 1:56) lines(1:10, pweib3(1:10, shp=weib90$shap[a], scl=weib90$scal[a], shift=weib90$shif[a]+weib90$offs[a]), col= mycolors28[a %% 28 + 1])
# abline(h=.9,lty=3)