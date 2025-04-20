
library(cir)
library(upndown)

# Dose-response function used in many examples
exampleMu = 5.6
exampleSig = 2
M = 10

# Logistic example used for most figures
exampleF = plogis(1:M, location = exampleMu, scale = exampleSig)
exampleF2 = plogis(1:(2*M), location = exampleMu, scale = exampleSig)

# Exponential "counter-example" :)
exp11F = pexp(1:M, rate=1/11)

outdir = '../../output'

# "Typical" plotting parameters

stdpar = par(mar = c(4,4,4,2), mgp = c(2.5, 0.6, 0), tck = -0.01, las = 1, cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.7)



