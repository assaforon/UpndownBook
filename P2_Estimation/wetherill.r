
wetherill <- function(x, y, rstart = 1, conf = 0.9, full = FALSE)
{
require(cir)
require(upndown)

revpts = reversals(y = y, x = x)
if (rstart > length(revpts) || is.na(revpts[1])) 
        return(NA)

# x at reversal points
revdoses = x[ revpts[rstart:length(revpts)] ]
# Point estimate is easy
pest = mean(revdoses)
if(is.null(conf)) return(pest)

# SE and CI
npairs = length(revdoses) %/% 2
pairmeans = ( revdoses[seq(1, 2*npairs-1, 2)] + revdoses[seq(2, 2*npairs, 2)] ) / 2
sterr = sd(pairmeans) / sqrt(npairs)

tmult = - qt((1-conf)/2, df = npairs-1)
confnames = paste(c('lower', 'upper'), round(100*conf), 'conf', sep='')

tmp = c(pest, pest - tmult*sterr, pest + tmult*sterr) 
names(tmp) = c("point", confnames)

if(!full) return(tmp)

return(list(revids = revpts, revdoses = revdoses, pairmeans = pairmeans, se = sterr, ests = tmp) )
}







	
