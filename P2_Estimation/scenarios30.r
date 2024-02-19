rm(list = ls() )
source('simulation_header.r')

#------------------------------ Constants / Utilities --------------------------#
set.seed(414155)
nsim = 1000
midtarg = 4.5
# Exclusion criteria and buffer factor
flat30 = 0.025
steep30 = 0.6
low30 = 0.05
high30 = 0.7
buff30 = 1.5

#### Note: the 'c' at the end of data table names indicated "candidate"
#          as some scenarios will be filtered out for being too steep/flat

#------------------------------ Weibull --------------------------#


weib30c = data.table(scal = c(runif(buff30*nsim/2, 0.5, 20), runif(buff30*nsim/2, 0.5, 6)),
					# shap = runif(buff30*nsim, .5, 10),
					shap = rgamma(buff30*nsim, shape = 0.8, scale = 5) + 0.5,
					offs = runif(buff30*nsim, -.5, .5) )
					
weib30c[ , shif := weibshift(shp=shap, scl=scal, targx = midtarg, targy = 0.3) ]

### F values at neighboring levels
weib30c[ , pp1 := pweib3(midtarg+0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib30c[ , pm1 := pweib3(midtarg-0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib30c[ , pm2 := pweib3(midtarg-1.5, shp=shap, scl=scal, shift = shif+offs) ]
					
# Targets: "exact..." (parametric)
weib30c[ , t30 := qweib3(0.3, shp=shap, scl=scal, shift = shif+offs) ]	
# ... and interpolated between closest levels (F-bar)				
weib30c[,row0 := 1:.N]
weib30c[ , t30bar := approx(c(pm1, pp1), midtarg + 0:1 - 0.5, xout = 0.3)$y, by = 'row0']


# Eligibility indicator
weib30c[ , elig := (pm1 > low30 & pp1 < high30 & 
						(pp1-pm1) %between% c(flat30, steep30) ) ]
weib30 = weib30c[weib30c$elig, ][1:nsim, ]
weib30[,row0 := 1:.N]

fwrite(weib30, file = file.path(outdir, 'scenarios_weib30.csv') )

mycolors36 <- c("antiquewhite4","aquamarine3","aquamarine4","azure4","black","blue","blueviolet","brown","burlywood4","chartreuse","chartreuse4","chocolate","coral2","coral4","cornflowerblue","darkgoldenrod","darkgoldenrod1","darkgreen","darkkhaki","darkmagenta","darkolivegreen","darkolivegreen2","darkorange2","darkred","darksalmon","darkslateblue","plum4","sienna2","slategray3","springgreen3","steelblue4","tan4","thistle4","tomato2","violet","violetred4")

stop('checkiout')

# mycolors28 <- c('red2','chartreuse2','gray35','cyan3','violetred3',mycolors36[-c(8,14,24,4,9,33,23,13,34,28,21,10,36)])
# mycolors28=mycolors28[c(9,1:4,10:13,16:28,15,5:8,14)]

# plot(c(1,12),0:1, type='n')
# for(a in 1:72) lines(1:12, pweib3(1:12, shp=weib30$shap[a], scl=weib30$scal[a], shift=weib30$shif[a]+weib30$offs[a]), col= mycolors36[a %% 28 + 1])
# abline(h=.9,lty=3)

#------------------------------ Logistic --------------------------#


logi30c = data.table(scal = 2^runif(buff30*nsim, -2, 3.),
					loc = runif(buff30*nsim, midtarg-.5, midtarg+.5) )
					
### F values at neighboring levels
logi30c[ , pp1 := plogis(midtarg+0.5, scale = scal, location =  loc) ]
logi30c[ , pm1 := plogis(midtarg-0.5, scale = scal, location =  loc) ]
logi30c[ , pm2 := plogis(midtarg-1.5, scale = scal, location =  loc) ]
					
# Targets: "exact..." (parametric)
logi30c[ , t30 := qlogis(0.5, scale = scal, location = loc) ]	
# ... and interpolated between closest levels (F-bar)				
logi30c[,row0 := 1:.N]
logi30c[ , t30bar := approx(c(pm1, pp1), midtarg + 0:1 - 0.5, xout = 0.5)$y, by = 'row0']


# Eligibility indicator
logi30c[ , elig := (pm1 > low30 & pp1 < high30 & 
						(pp1-pm1) %between% c(flat30, steep30) ) ]

logi30 = logi30c[logi30c$elig, ][1:nsim, ]
logi30[,row0 := 1:.N]

fwrite(logi30, file = file.path(outdir, 'scenarios_logi30.csv') )

