source('basics_header.r')
library(data.table)
library(ggplot2)
theme_set(theme_bw(20))

groupsizes = 1:6

groups = CJ(K = groupsizes, l = groupsizes-1, u = groupsizes)
# Selecting only "legal" (K,l,u) trios
groups = groups[l < u & u <= K & K > 1, ]

groups[ , fstar := g2targ(K, l, u), by = .I  ]
groups[ , trio := paste('(', K, ',', l,',', u,')', sep=''), by = .I ] 

p1 <- ggplot(groups, aes(x=fstar, y = K + (u - l - K/2 - .5)/10) ) + geom_text(aes(label = trio), size = 6)  +
		scale_x_continuous(breaks = seq(.1,.9,.1) ) + scale_y_continuous(breaks = 1:6 ) +
		labs(x = "Balance Point's Response Rate (F*)", y = 'Group Size (K)') + 
		theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
		axis.ticks.length.y  = unit(0., "cm") )

ggsave(p1, file = file.path(outdir, 'gudTargets.pdf'), width = 12, height = 7) 
