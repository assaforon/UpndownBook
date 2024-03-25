#### Doing all the estimates after n=30 instead of 60
cat(base::date(), '\n')
rm(list=ls())
outdir = '../../output'

load(file.path(outdir, 'grandsim30w.RData'))	
load(file.path(outdir, 'grandsim30w_gud.RData'))	
source('simulation_header.r')

nsims = 500


library(magrittr)

hestkw30minmid = estbatch(kw30minmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=krow, curvy = TRUE, desargs=k30list)
	
hestkw30midmid = estbatch(kw30midmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=krow, curvy = TRUE, desargs=k30list)

hestkw30minhi = estbatch(kw30minhi, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=krow, curvy = TRUE, desargs=k30list)
hestkw30minlo = estbatch(kw30minlo, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=krow, curvy = TRUE, desargs=k30list)
	
cat('k standard\n')	
					
########### BCD standard	

hestbw30minmid = estbatch(bw30minmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, 
					 desfun=bcd, curvy = TRUE, desargs=b30list)
		
hestbw30midmid = estbatch(bw30midmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3,  
				desfun=bcd, curvy = TRUE, desargs=b30list)

hestbw30minhi = estbatch(bw30minhi, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					 desfun=bcd, curvy = TRUE, desargs=b30list)

hestbw30minlo = estbatch(bw30minlo, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					 desfun=bcd, curvy = TRUE, desargs=b30list)

cat('b standard\n')						
									
#### GUD

hestg2w30minmid = estbatch(g2w30minmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, 
					bpt=ktarg30, desfun=groupUD, curvy = TRUE, desargs=g2list)

hestg2w30midmid = estbatch(g2w30midmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, bpt=ktarg30, 
				desfun=groupUD, curvy = TRUE, desargs=g2list)

hestg2w30minhi = estbatch(g2w30minhi, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					bpt=ktarg30, desfun=groupUD, curvy = TRUE, desargs=g2list)

hestg2w30minlo = estbatch(g2w30minlo, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]-2, target=.3, 
					bpt=ktarg30, desfun=groupUD, curvy = TRUE, desargs=g2list)
		
cat('GUD 2\n')						
				
hestg3w30minmid = estbatch(g3w30minmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3, 
					 desfun=groupUD, curvy = TRUE, desargs=g3list, bpt = gtarg302)

hestg3w30midmid = estbatch(g3w30midmid, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim], target=.3,  
				desfun=groupUD, curvy = TRUE, desargs=g3list, bpt = gtarg302)

hestg3w30minhi = estbatch(g3w30minhi, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]+2, target=.3, 
					 desfun=groupUD, curvy = TRUE, desargs=g3list, bpt = gtarg302)

hestg3w30minlo = estbatch(g3w30minlo, n=30, nsim=nsims, truth=weib30parm$t30[1:nsim]-2, target=.3, 	 desfun=groupUD, curvy = TRUE, desargs=g3list, bpt = gtarg302)

save.image(file.path(outdir, 'halfsim30w.RData'))					
cat('GUD 3 and done.\n')					
									



		

		
					





