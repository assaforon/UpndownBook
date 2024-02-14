cat(base::date(), '\n')
source('simulation_header.r')


### Constants

pointnames = c("dm48", "all1", "all3", "rev1", "dyna", "cir", "ir" )

intnames = c("cir", "ir", "dyna", "all3", "rev1", "cboot")


#-------------------------- Atomic Utilities

### Point estimate performance
pmetrix <- function(simout, estnames = pointnames, bigerr = 0.95) 
{
	require(data.table)
	calcnames = intersect(estnames, names(simout) )

simout[ ,apply(.SD, 2, trio, ref=true, p=bigerr), .SDcol = calcnames ]
}

### Interval performance (coverage and width)

imetrix <- function(simout)
{
	require(data.table)
	simdat = copy(simout)
	
dout = simdat[ , list(all3 = mean(all3l<=true & all3u>=true, na.rm=TRUE),
	rev1 = mean(rev1l<=true & rev1u>=true, na.rm=TRUE),
	dyna = mean(dynal<=true & dynau>=true, na.rm=TRUE),
	cir = mean(cirl<=true & ciru>=true, na.rm=TRUE),
	ir = mean(irl<=true & iru>=true, na.rm=TRUE),
	cirboot = mean(cbootl<=true & cbootu>=true, na.rm=TRUE)
	) ]

	simdat[ , cirfin := (is.finite(ciru) & is.finite(cirl) ) ]
	simdat[ , irfin :=  (is.finite(iru) & is.finite(irl) ) ]

rbind(dout, simdat[ , list(all3 = mean(all3u - all3l, na.rm=TRUE),
	rev1 = mean(rev1u - rev1l, na.rm=TRUE),
	dyna = mean(dynau - dynal, na.rm=TRUE),
	cir = mean(ciru[cirfin] - cirl[cirfin]),
	ir = mean(iru[irfin] - irl[irfin]),
	cfinite = mean(cirfin) , ifinite = mean(irfin),
#	liao = ifelse(addLiao, mean(liaou - liaol, na.rm=TRUE), NA),
	cirboot = mean(cbootu - cbootl, na.rm=TRUE) ) ], fill = TRUE )
	
}
	