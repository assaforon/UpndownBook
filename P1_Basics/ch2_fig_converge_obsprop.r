source('basics_header.r')

### Additional constants/utilities
bpi = pivec(exampleF, bcdmat, target = btarg)

# Cumulative allocation frequency of dose m
cumulp = function(x, m) cumsum(x==m)/(1:m) 


p4 = apply(longdat$doses, 2, cumulp, m=4)
p1 = apply(longdat$doses, 2, cumulp, m=1)

