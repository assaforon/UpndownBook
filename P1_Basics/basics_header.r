
library(cir)
library(upndown)

# Dose-response function used in many examples
exampleMu = 5.6
exampleSig = 2
M = 10

# Logistic example used for most figures
exampleF = plogis(1:M, location = exampleMu, scale = exampleSig)
# Exponential "counter-example" :)
exp11F = pexp(1:M, rate=1/11)

outdir = '../../output'



