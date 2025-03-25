source('basics_header.r')

### Additional constants/utilities
btarg = 0.3
n = 5000
nsim = 200
x1 = 2

set.seed(1994)

### Running a long-n simulation to get a distribution of first-visit times

# BCD - logistic
blongdat1 = dfsim(n, Fvals = exampleF, starting = x1, design = bcd, 
			desArgs = list(coin = btarg/(1-btarg), lowTarget=TRUE), ensemble = nsim)

saveRDS(blongdat1, file = file.path(outdir, 'ch3_sim_firstvisit_distribution_bcd1.rds') )

# BCD - exponential
blongdat2 = dfsim(n, Fvals = exp11F, starting = x1, design = bcd, 
			desArgs = list(coin = btarg/(1-btarg), lowTarget=TRUE), ensemble = nsim)

saveRDS(blongdat2, file = file.path(outdir, 'ch3_sim_firstvisit_distribution_bcd2.rds') )

# Classical (which is the default of dfsim(), so no design arguments needed)
clongdat = dfsim(n, Fvals = exampleF, starting = x1, ensemble = nsim)

saveRDS(clongdat, file = file.path(outdir, 'ch3_sim_firstvisit_distribution_classical.rds') )
