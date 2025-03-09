source('basics_header.r')


cpi = pivec(exampleF, classicmat)

eye = diag(1, M, M) # identity matrix I
Z = solve(eye - classicmat(exampleF) + t(matrix(rep(cpi, M), nrow = M)) )

startdose = 2

firstimes = ( diag(Z) - Z[startdose, ] ) / cpi
firstsq =  ( 2 * diag(Z) - cpi ) / (cpi^2)
firstsd = sqrt(firstsq - firstimes^2)

# Kemeny and Snell
# f_ij=number of steps before reaching d_j, starting at d_i.
#w_ij= E_i[f_j^2]  = -(1/pi_j) + (2z_jj) / (pi_jj)^2