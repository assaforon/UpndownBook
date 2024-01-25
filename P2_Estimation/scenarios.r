source('simulation_header.r')

#------------------------------ Constants / Utilities --------------------------#
set.seed(414155)
nsim = 1000
midtarg = 6.5
# Exclusion criteria and buffer factor
flat90 = 0.03
low90 = 0.5
buff90 = 1.5


#------------------------------ 90th percentile --------------------------#

weib90c = data.table(scal = runif(buff90*nsim, 2, 10),
					# shap = runif(buff90*nsim, .5, 10),
					shap = rgamma(buff90*nsim, shape = 2., scale = 1.5) + 0.8,
					offs = runif(buff90*nsim, -.5, .5) )
					
weib90c[ , shif := weibshift(shp=shap, scl=scal, targx = midtarg, targy = 0.9) ]

### F values at neighboring levels
weib90c[ , pp1 := pweib3(midtarg+0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib90c[ , pm1 := pweib3(midtarg-0.5, shp=shap, scl=scal, shift = shif+offs) ]
weib90c[ , pm2 := pweib3(midtarg-1.5, shp=shap, scl=scal, shift = shif+offs) ]
					
# Targets: "exact..." (parametric)
weib90c[ , t90 := qweib3(0.9, shp=shap, scl=scal, shift = shif+offs) ]	
# ... and interpolated between closest levels (F-bar)				
weib90c[,row0 := 1:.N]
weib90c[ , t90bar := approx(c(pm1, pp1), midtarg + 0:1 - 0.5, xout = 0.9)$y, by = 'row0']


# Eligibility indicator
weib90c[ , elig := (pm1 > low90 & pp1-pm1 > flat90) ]
weib90 = weib90c[weib90c$elig, ][1:nsim, ]

fwrite(weib90, file = file.path(outdir, 'scenarios_weib90.csv') )

# mycolors28 <- c('red2','chartreuse2','gray35','cyan3','violetred3',mycolors36[-c(8,14,24,4,9,33,23,13,34,28,21,10,36)])
# mycolors28=mycolors28[c(9,1:4,10:13,16:28,15,5:8,14)]

# plot(c(1,12),0:1, type='n')
# for(a in 1:56) lines(1:12, pweib3(1:12, shp=weib90$shap[a], scl=weib90$scal[a], shift=weib90$shif[a]+weib90$offs[a]), col= mycolors28[a %% 28 + 1])
# abline(h=.9,lty=3)