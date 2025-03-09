source('basics_header.r')

### Additional constants/utilities
btarg = 0.3
n = 2000
nsim = 500

### Running a long-n simulation

longdat = dfsim(n, Fvals = exampleF, starting = 1, design = bcd, 
			desArgs = list(coin = btarg/(1-btarg), lowTarget=TRUE), ensemble = nsim)

saveRDS(longdat, file = file.path(outdir, 'ch2_sim_converge_obsprop.rds') )
