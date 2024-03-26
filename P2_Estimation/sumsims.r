#-----------------------------------------------------------
#---------  Summary and viz functions for simulation output
#-----------------------------------------------------------


cat(base::date(), '\n')
source('simulation_header.r')


### Constants

pointnames = c("dm48", "rev1", "all1", "all3", "dyna", "cir", "ir", "dual")
pointnice = c('Dixon-Mood', 'Revs. (Wetherill)', paste('Avg. from' ,c('R1','R3') ),  
			'Dynamic Avg.', 'CIR', 'IR', "Combination" )

intnames = c("rev1", "all3", "dyna", "cir", "ir", 'dual', "cirboot")
intnice = c('Revs. (Wetherill)', 'Avg. from R3', 'Dynamic Avg.', 'CIR', 'IR' ,   "Combination", 'CIR Boot.')


#-------------------------- Atomic Utilities

### Point estimate performance
pmetrix <- function(simout, estnames = pointnames, combine = TRUE, 
					bigerr = 0.9, apples2apples = TRUE,...) 
{
require(data.table)
simdat = copy(simout)

# combining CIR and dynamean
if(combine) simdat[ , dual := ifelse(is.na(cir), dyna, (cir+dyna)/2 ) ]
# Potentially excluding rows w/o CIR estimate
if(apples2apples) simdat = simdat[!is.na(cir), ]

calcnames = intersect(estnames, names(simdat) )

simdat[ ,apply(.SD, 2, trio, ref=true, p=bigerr), .SDcol = calcnames ]
}

### Interval performance (coverage and width)

imetrix <- function(simout, combine = TRUE, apples2apples = TRUE,  ...)
{
require(data.table)
simdat = copy(simout)

# combining CIR and dynamean
if(combine)
{
	simdat[ , duall := ifelse(is.na(cirl), dynal, (cirl+dynal)/2 ) ]
	simdat[ , dualu := ifelse(is.na(ciru), dynau, (ciru+dynau)/2 ) ]
} 
# Potentially excluding rows w/o CIR estimate
if(apples2apples) simdat = simdat[!is.na(cir), ]
	
dout = simdat[ , list(all3 = mean(all3l<=true & all3u>=true, na.rm=TRUE),
	rev1 = mean(rev1l<=true & rev1u>=true, na.rm=TRUE),
	dyna = mean(dynal<=true & dynau>=true, na.rm=TRUE),
	cir = mean(cirl<=true & ciru>=true, na.rm=TRUE),
	dual = ifelse(combine, mean(duall<=true & dualu>=true, na.rm=TRUE), NA),
	ir = mean(irl<=true & iru>=true, na.rm=TRUE),
	cirboot = ifelse(exists('cbootl'), mean(cbootl<=true & cbootu>=true, na.rm=TRUE), NA) )
	 ]

	simdat[ , cirfin := (is.finite(ciru) & is.finite(cirl) ) ]
	simdat[ , irfin :=  (is.finite(iru) & is.finite(irl) ) ]
#	simdat[ , cbfin :=  ifelse(exists('cbootl'), 
#						(is.finite(cbootu) & is.finite(cbootl) ), NA) ]


rbind(dout, simdat[ , list(all3 = mean(all3u - all3l, na.rm=TRUE),
	rev1 = mean(rev1u - rev1l, na.rm=TRUE),
	dyna = mean(dynau - dynal, na.rm=TRUE),
	cir = mean(ciru[cirfin] - cirl[cirfin]),
	dual = ifelse(combine, mean(dualu - duall, na.rm=TRUE), NA),
	ir = mean(iru[irfin] - irl[irfin])
#	cfinite = mean(cirfin) 
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

combo <- function(simlist, atomfun = pmetrix, finites = TRUE,...) 
{
	require(data.table)
	
	metrics = lapply(simlist, function(x) atomfun(get(x), ...) )
	
	dout = data.table()
	for(a in seq_along(simlist) ) 
	{
#	cat(simlist[a],'\n')
		dout = rbind(dout, stackem(metrics[[a]], 
			framework = gsub('est', '', simlist[a]), ... ) )
			}
	if(finites)
	{
		tmp = sapply(simlist, function(x) mean( is.finite(get(x)$cir) ) )
		dout = rbind(dout, data.table(Metric = 'Finite CIR', 
		Framework = gsub('est', '', names(tmp) ), Value = tmp), fill = TRUE)
	}
	dout
}

#--------------------------  Plotting

sideside <- function(omnibus, metric = 'RMSE', fsize = 15, jwid = 0.1, yref = NULL,
			zoom = c(0, NA), expansion = c(0, 0.02), psize = 4, rotlab = TRUE,
			innames = pointnames, outnames = pointnice, desvar = TRUE, multip = 1,
			colkey = c('grey65', 'black'), titl = '', addmean = TRUE, ytext =  NULL)
{
require(plyr)
require(forcats)
require(data.table)
require(ggplot2)
theme_set(theme_bw(fsize))

pdat = omnibus[Metric == metric & estimate %in% innames, ]
pdat[ , Estimate := fct_relevel(mapvalues(estimate, innames, outnames, 
						warn_missing = FALSE), outnames) ]
						
if(addmean) pmeans = pdat[, list(Value = mean(multip*Value)), keyby = 'Estimate' ]

pout = ggplot(pdat, aes(Estimate, multip*Value))  +
                scale_y_continuous(expand = expansion, limits = zoom )	+ 
				labs(y = ifelse(is.null(ytext), paste(metric, '(spacing units)'), ytext), x='', title = titl)

if(!is.null(yref)) pout = pout + geom_hline(yintercept = yref)

if(addmean) pout = pout + geom_point(data=pmeans, aes(Estimate, Value), pch='-', size = psize*8 )

if(desvar) {
	pout = pout + geom_jitter(width = jwid, size=psize, 
						aes(Estimate, multip*Value, color = Design) ) + 
						scale_color_manual(values = colkey) } else {
		pout =	pout + geom_jitter(width = jwid, size=psize) }
		
if(rotlab) pout = pout + theme(axis.text.x = element_text(angle = 45, hjust=1) )

pout
}


				



	