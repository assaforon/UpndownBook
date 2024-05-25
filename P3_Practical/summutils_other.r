
### Extracting metrics from a list of ensembles

distill <- function(simlist, yfun = function(x) sum(x$goodest, na.rm = TRUE)/nrow(x), xfun = function(x) mean(x$responses), cuechar = 'w')
{
require(stringr)
require(data.table)

estlist = sort(grep('rest', simlist, value = TRUE))
runlist = sort(setdiff(simlist, estlist))

dout = data.table(Framework = runlist)
dout[ , des := tstrsplit(Framework, split = cuechar)[[1]] ]
dout[ , des := str_trim(paste(des, gsub('[a-z]', '', str_sub(Framework, -5)) ) ) ]

dout$y = sapply(estlist, function(x) yfun(get(x)) )  
dout$x = sapply(runlist, function(x) xfun(get(x)) )  

return(dout)
}

### Extracting a ensemble data from name, and adding names as variable

dextract <- function(dname, cuechar1 = 'w', cuechar2 = 0)
{
require(stringr)
require(data.table)

dout = get(dname)
dout[ , Framework := gsub('rest', '', dname) ]
dout[ , des := tstrsplit(Framework, split = cuechar1)[[1]] ]
dout[ , des := str_trim(paste(des, gsub('[a-z]', '', str_sub(Framework, -5)) ) ) ]
dout[ , sett := tstrsplit(Framework, split = cuechar2)[[2]] ]
dout[ , sett := gsub('[0-9]', '', sett) ]
dout[ , sett := gsub('[_]', '', sett) ]

dout
}



# , xref, yref = 0.3