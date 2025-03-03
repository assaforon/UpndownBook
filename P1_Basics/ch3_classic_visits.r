source('basics_header.r')


cpi = pivec(exampleF, classicmat)

eye = diag(1, M, M) # identity matrix I
Z = solve(eye - classicmat(exampleF) + t(matrix(rep(cpi, M), nrow = M)) )

startdose = 2

firstimes = ( diag(Z) - Z[startdose, ] ) / cpi