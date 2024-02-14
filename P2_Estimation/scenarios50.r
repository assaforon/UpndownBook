rm(list = ls() )
source('simulation_header.r')

#------------------------------ Constants / Utilities --------------------------#
set.seed(414155)
nsim = 1000
midtarg = 6.5
# Exclusion criteria and buffer factor
flat50 = 0.025
steep50 = 0.6
low50 = 0.1
high50 = 0.9
buff50 = 1.5

#### Note: the 'c' at the end of data table names indicated "candidate"
#          as some scenarios will be filtered out for being too steep/flat

#------------------------------ Weibull --------------------------#


weib50c = data.table(scal = runif(buff50*nsim, 1, 15),
					# shap = runif(buff50*nsim, .5, 10),
					shap = rgamma(buff50*nsim, shape = 0.7, scale = 6) + 0.5,
					offs = runif(buff50*nsim, -.5, .5) )
					
weib50c[ , shif := weibshift(shp=shap, scl=scal, targx = midtarg, targy = 0.5) ]

### F values at neighboring levels
weib50c[ , pp1 := pweib3(midtarg+0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib50c[ , pm1 := pweib3(midtarg-0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib50c[ , pm2 := pweib3(midtarg-1.5, shp=shap, scl=scal, shift = shif+offs) ]
					
# Targets: "exact..." (parametric)
weib50c[ , t50 := qweib3(0.5, shp=shap, scl=scal, shift = shif+offs) ]	
# ... and interpolated between closest levels (F-bar)				
weib50c[,row0 := 1:.N]
weib50c[ , t50bar := approx(c(pm1, pp1), midtarg + 0:1 - 0.5, xout = 0.5)$y, by = 'row0']


# Eligibility indicator
weib50c[ , elig := (pm1 > low50 & pp1 < high50 & 
						(pp1-pm1) %between% c(flat50, steep50) ) ]
weib50 = weib50c[weib50c$elig, ][1:nsim, ]
weib50[,row0 := 1:.N]

fwrite(weib50, file = file.path(outdir, 'scenarios_weib50.csv') )

mycolors36 <- c("antiquewhite4","aquamarine3","aquamarine4","azure4","black","blue","blueviolet","brown","burlywood4","chartreuse","chartreuse4","chocolate","coral2","coral4","cornflowerblue","darkgoldenrod","darkgoldenrod1","darkgreen","darkkhaki","darkmagenta","darkolivegreen","darkolivegreen2","darkorange2","darkred","darksalmon","darkslateblue","plum4","sienna2","slategray3","springgreen3","steelblue4","tan4","thistle4","tomato2","violet","violetred4")

#stop('checkiout')

# mycolors28 <- c('red2','chartreuse2','gray35','cyan3','violetred3',mycolors36[-c(8,14,24,4,9,33,23,13,34,28,21,10,36)])
# mycolors28=mycolors28[c(9,1:4,10:13,16:28,15,5:8,14)]

# plot(c(1,12),0:1, type='n')
# for(a in 1:72) lines(1:12, pweib3(1:12, shp=weib50$shap[a], scl=weib50$scal[a], shift=weib50$shif[a]+weib50$offs[a]), col= mycolors36[a %% 28 + 1])
# abline(h=.9,lty=3)

#------------------------------ Logistic --------------------------#


logi50c = data.table(scal = 2^runif(buff50*nsim, -2, 3.),
					loc = runif(buff50*nsim, midtarg-.5, midtarg+.5) )
					
### F values at neighboring levels
logi50c[ , pp1 := plogis(midtarg+0.5, scale = scal, location =  loc) ]
logi50c[ , pm1 := plogis(midtarg-0.5, scale = scal, location =  loc) ]
logi50c[ , pm2 := plogis(midtarg-1.5, scale = scal, location =  loc) ]
					
# Targets: "exact..." (parametric)
logi50c[ , t50 := qlogis(0.5, scale = scal, location = loc) ]	
# ... and interpolated between closest levels (F-bar)				
logi50c[,row0 := 1:.N]
logi50c[ , t50bar := approx(c(pm1, pp1), midtarg + 0:1 - 0.5, xout = 0.5)$y, by = 'row0']


# Eligibility indicator
logi50c[ , elig := (pm1 > low50 & pp1 < high50 & 
						(pp1-pm1) %between% c(flat50, steep50) ) ]

logi50 = logi50c[logi50c$elig, ][1:nsim, ]
logi50[,row0 := 1:.N]

fwrite(logi50, file = file.path(outdir, 'scenarios_logi50.csv') )

