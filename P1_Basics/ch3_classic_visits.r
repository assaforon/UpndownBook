source('basics_header.r')

startdose = 2

cpi = pivec(exampleF, classicmat)

eye = diag(1, M, M) # identity matrix I
Z = solve(eye - classicmat(exampleF) + t(matrix(rep(cpi, M), nrow = M)) )
Zdg = diag(diag(Z)) # version of Z with 0 off-diags
# K & S book matrices
dee = diag(1/cpi, M, M) # D
E = matrix(1, M, M) 

initvec = rep(0, M)
initvec[startdose] = 1

# Mean first-passage time per K & S
Mij =  (eye - Z + E %*% Zdg ) %*% dee
fij = initvec %*% Mij
# Square of that
Wij = Mij %*% (2 * Zdg %*% dee - eye) + 2* (Z %*% Mij - E %*% diag(diag(Z %*% Mij)) )
f2ij = initvec %*% Wij

fsd = sqrt(f2ij - fij^2)

## Distribution of visits

classic5k = readRDS(file.path(outdir, 'ch3_sim_firstvisit_distribution_classical.rds') ) 

firstvisit = function(x, m) min(which(x == m))

firsts = matrix(NA, nrow = M, ncol = ncol(classic5k$doses) )

for(m in 1:M) firsts[m,] = apply(classic5k$doses[-1, ], 2, firstvisit, m=m)



#firstimes = ( diag(Z) - Z[startdose, ] ) / cpi
#firstsq =  ( 2 * diag(Z) - cpi ) / (cpi^2)
#firstsd = sqrt(firstsq - firstimes^2)

# Kemeny and Snell
# f_ij=number of steps before reaching d_j, starting at d_i.
#w_ij= E_i[f_j^2]  = -(1/pi_j) + (2z_jj) / (pi_jj)^2