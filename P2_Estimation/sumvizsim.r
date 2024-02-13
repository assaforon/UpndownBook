cat(base::date(), '\n')
source('simulation_header.r')

### Utilities



metrix <- function(simout, bigerr = 0.95)

### Point estimate performance
tmp=list(metrics=ests[ ,apply(.SD, 2, trio, ref=true, p=bigerr), .SDcol=names(ests)[2:8] ],
	cirmissed=mean(is.na(ests$cir)), n = n, ensemble = nsim, 
	target = target, startpt = simdat$doses[1,1], targloc = mean(truth))

### CI coverage

tmp$coverage = ests[ , list(all3 = mean(all3l<=true & all3u>=true, na.rm=TRUE),
	rev1 = mean(rev1l<=true & rev1u>=true, na.rm=TRUE),
	dyna = mean(dynal<=true & dynau>=true, na.rm=TRUE),
	cir = mean(cirl<=true & ciru>=true, na.rm=TRUE),
	ir = mean(irl<=true & iru>=true, na.rm=TRUE),
	cirboot = mean(cbootl<=true & cbootu>=true, na.rm=TRUE)
	) ]

### CI width
	ests[ , cirfin := (is.finite(ciru) & is.finite(cirl) ) ]
	ests[ , irfin :=  (is.finite(iru) & is.finite(irl) ) ]
	
tmp$widths = ests[ , list(all3 = mean(all3u - all3l, na.rm=TRUE),
	rev1 = mean(rev1u - rev1l, na.rm=TRUE),
	dyna = mean(dynau - dynal, na.rm=TRUE),
	cir = mean(ciru[cirfin] - cirl[cirfin]),
	ir = mean(iru[irfin] - irl[irfin]),
	cfinite = mean(cirfin) , ifinite = mean(irfin),
#	liao = ifelse(addLiao, mean(liaou - liaol, na.rm=TRUE), NA),
	cirboot = mean(cbootu - cbootl, na.rm=TRUE) ) ]