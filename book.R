### ================================================================= ###
## GAUSSIAN MARKOV RANDOM FIELDS ===================================== ##
## Theory and Applications =========================================== ##
## Haavard Rue & Leonhard Held ======================================= ##
### ================================================================= ###

### Chapter 2: Theory of Gaussian Markov random fields ==================
## 2.3 Simulation from a GRMF ===========================================

# Algorithm 2.1 =========================================================
##  Solving \mathbf{Ax} = \mathbf{b} where \mathbf{A} > 0 ============ ##

da <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 5), sigma = diag(5))
( A <- cov(da) )

( b <- as.matrix( rnorm(5) ) )

## 1: Compute the Cholesky factorization, \mathbf{A} = \mathbf{LL}^{T} ##
( esky <- chol(A) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
L <- t(esky) ; tL <- esky

round( L %*% tL, 15 ) == round( A, 15 )

## 2: Solve \mathbf{Lv} = \mathbf{b} ========== (Forward substitution) ##
diag(L) == diag(tL)
diag(L) != 0

( v <- numeric(5) )

v[1] <- b[1] / L[1, 1]

for (i in 2:5)
  v[i] <- ( 1 / L[i, i] ) * ( b[i] - sum( L[i, 1:(i-1)] * v[1:(i-1)] ) )
v

## 3: Solve \mathbf{L}^{T}\mathbf{x} = \mathbf{v}  (Back substitution) ##
( x <- numeric(5) )

x[5] <- v[5] / L[5, 5]

for (i in 4:1)
  x[i] <- ( 1 / L[i, i] ) * ( v[i] - sum( L[(i+1):5, i] * x[(i+1):5] ) )
## 4: Return \mathbf{x} ============================================== ##
x

cbind( A %*% x, b )
round( A %*% x, 14 ) == round( b, 14 )
### ================================================================= ###