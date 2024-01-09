rm(list=ls())
library(upndown)

# For use on your own machine, change 'outdir' to where you want the figures saved
outdir = ('../../output')



## Boundary Effect
library(tis) # for barplot2(), shifted barplots

m=8
sfac=250
axfac=120
labfac=100
barcol='blue'
barwid=11
hashdens=15 # boundary density
hashcol='gray50'
hashwid=0.25 # boundary width

ymax=38
targ=4.8
sig=2
figF=plogis(1:m,location=targ,scale=sig)
exc=8

pdf(file.path(outdir,'averaging1.pdf'),width = 8, height = 14)
par(mfrow=c(4,1), tcl=-0.4, mar=c(1,6.5,1,1), mgp=c(4,1.,0), las=1)

pis=list()
means=rep(NA,4)
for (a in 0:3)
{
	pis[[a+1]]=pivec(cdf = plogis((1:m)-a,location=targ,scale=sig), matfun = classicmat)
	means[a+1]=weighted.mean((1:m)-a,w=pis[[a+1]])

	plot((1:m)-a,100*pis[[a+1]],type='h',xlim=c(-2.3,m+0.3),
	xlab=ifelse(a<3,'','Dose'),ylab='Percent of Patients',ylim=c(0,ymax),xaxt='n',
	cex.lab=sfac/labfac,cex.axis=sfac/axfac,lwd=barwid,col=barcol)
	abline(v=targ,lwd=3, lty = 2)
	abline(v=means[a+1],lwd=3)

# Boundaries
	par(lwd=3,yaxt='n')
	barplot2(ymax,x.offset=m-a+0.2,width=hashwid,density=hashdens,add=TRUE,border=NA,angle=30,col=hashcol)
	barplot2(ymax,x.offset=m-a+0.2,width=hashwid,density=hashdens,add=TRUE,border=NA,angle=-30,col=hashcol)
	if(a<4) barplot2(ymax,x.offset=0.6-a,width=hashwid,density=hashdens,add=TRUE,border=NA,angle=30,col=hashcol)
	if(a<4)barplot2(ymax,x.offset=0.6-a,width=hashwid,density=hashdens,add=TRUE,border=NA,angle=-30,col=hashcol)
	par(lwd=1,yaxt='s')

}

 dev.off()
