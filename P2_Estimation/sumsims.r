#-----------------------------------------------------------
#---------  Summary and viz functions for simulation output
#-----------------------------------------------------------


cat(base::date(), '\n')
source('simulation_header.r')


### Constants

pointnames = c("dm48", "all1", "all3", "rev1", "dyna", "cir", "ir" )
pointnice = c('Dixon-Mood', paste('Avg. from' ,c('R1','R3') ), 'Revs. (Wetherill)', 
			'Dynamic Avg.', 'CIR', 'IR' )

intnames = c("cir", "ir", "dyna", "all3", "rev1", "cboot")
intnice = c('CIR', 'IR', 'Dynamic Avg.', 'Avg. from R3', 
			'Revs. (Wetherill)', 'CIR Boot.' )


#-------------------------- Atomic Utilities

### Point estimate performance
pmetrix <- function(simout, estnames = pointnames, bigerr = 0.95, ...) 
{
	require(data.table)
	calcnames = intersect(estnames, names(simout) )

simout[ ,apply(.SD, 2, trio, ref=true, p=bigerr), .SDcol = calcnames ]
}

### Interval performance (coverage and width)

imetrix <- function(simout, ...)
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
	simdat[ , cbfin :=  (is.finite(cbootu) & is.finite(cbootl) ) ]


rbind(dout, simdat[ , list(all3 = mean(all3u - all3l, na.rm=TRUE),
	rev1 = mean(rev1u - rev1l, na.rm=TRUE),
	dyna = mean(dynau - dynal, na.rm=TRUE),
	cir = mean(ciru[cirfin] - cirl[cirfin]),
	ir = mean(iru[irfin] - irl[irfin]),
	cfinite = mean(cirfin) 
#	liao = ifelse(addLiao, mean(liaou - liaol, na.rm=TRUE), NA),
#	cirboot = mean(cbootu[cbfin] - cbootl[cbfin], na.rm=TRUE) 
	) ], fill = TRUE )
}

### Stackem up!

stackem <- function(metr, outnames = NULL, framework, finites = TRUE, ...)
{
	require(data.table)
	if(!is.null(outnames)) rownames(metr) = outnames
	
	tmp = data.table(metr)
	tmp[ , Metric := rownames(metr) ]
	
	dout = melt(tmp, id.var = 'Metric', variable.name = 'estimate', 
			value.name = 'Value', na.rm = finites, variable.factor = FALSE)
	dout[ , Framework := framework ]
	dout
}

#-------------------------- Combining Atomic Utilities

combo <- function(simlist, atomfun = pmetrix, ...) 
{
	require(data.table)
	
	metrics = lapply(simlist, function(x) atomfun(get(x), ...) )
	
	dout = data.table()
	for(a in seq_along(simlist) ) 
	{
#	cat(simlist[a],'\n')
		dout = rbind(dout, stackem(metrics[[a]], 
			framework = gsub('^est', '', simlist[a]), ... ) )
			}
	dout
}

#--------------------------  Plotting

sideside <- function(omnibus, metric = 'RMSE', fsize = 15, jwid = 0.1, yref = NULL,
			zoom = c(0, NA), expansion = c(0, 0.01), psize = 3,
			innames = pointnames, outnames = pointnice, colvar = FALSE,
			colkey = c('grey65', 'black'), titl = '' )
{
require(plyr)
require(data.table)
require(ggplot2)
theme_set(theme_bw(fsize))

pdat = omnibus[Metric == metric & estimate %in% innames, ]
pdat[ , Estimate := mapvalues(estimate, innames, outnames) ]

pout = ggplot(pdat, aes(Estimate, Value))  +
                scale_y_continuous(expand = expansion, limits = zoom )	+ 
				labs(y = metric, x='', title = titl)

if(!is.null(yref)) pout = pout + geom_hline(yintercept = yref)

if(colvar) {
	pout = pout + geom_jitter(width = jwid, size=psize, 
						aes(Estimate, Value, color = Design) ) + 
						scale_color_manual(values = colkey) } else {
		pout =	pout + geom_jitter(width = jwid, size=psize) }

pout
}


				



	