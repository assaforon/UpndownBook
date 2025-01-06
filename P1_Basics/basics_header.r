
library(cir)
library(upndown)

# Dose-response function used in many examples
exampleMu = 5.6
exampleSig = 2
exampleF = plogis(1:10, location = exampleMu, scale = exampleSig)

outdir = '../../output'

