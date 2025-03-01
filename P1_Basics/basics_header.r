
library(cir)
library(upndown)

# Dose-response function used in many examples
exampleMu = 5.6
exampleSig = 2
M = 10
exampleF = plogis(1:M, location = exampleMu, scale = exampleSig)

outdir = '../../output'

