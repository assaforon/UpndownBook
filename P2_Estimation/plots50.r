cat(base::date(), '\n')
rm(list=ls())
library(data.table)

outdir = '../../output'

wid = 12
hgt = 7

load(file.path(outdir, 'plotting.RData'))	

library(ggplot2)
theme_set(theme_bw(16))

ggsave(point50r + labs(y = "RMSE (spacing units)", title = ''),
			file = file.path(outdir, 'sim_rmse50.pdf'), width=wid, height=hgt)

ggsave(point50b + labs(y = "Bias (spacing units)", title = ''),
			file = file.path(outdir, 'sim_bias50.pdf'), width=wid, height=hgt)

ggsave(int50c + labs(y = "Interval Coverage (%)", title = ''),
			file = file.path(outdir, 'sim_cover50.pdf'), width=wid, height=hgt)

ggsave(int50w + labs(y = "Average Interval Width (spacing units)", title = ''),
			file = file.path(outdir, 'sim_width50.pdf'), width=wid, height=hgt)
			
			
cat(base::date(), '\n')





