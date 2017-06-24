### ================================================================= ###
## GAUSSIAN MARKOV RANDOM FIELDS ===================================== ##
## Theory and Applications =========================================== ##
## Haavard Rue & Leonhard Held ======================================= ##
### ================================================================= ###

### Chapter 2: Theory of Gaussian Markov random fields ==================
## 2.3 Simulation from a GRMF ===========================================
## 2.3.1 Some basic numerical linear algebra ============================

# Algorithm 2.1 =========================================================
##  Solving \mathbf{Ax} = \mathbf{b} where \mathbf{A} > 0 ============ ##

da <- mvtnorm::rmvnorm(n = 10, mean = rep(0, 5), sigma = diag(5))
( A <- cov(da) ) # A is SPB
## or:
# df <- data.frame( a = gl(3, 4), b = gl(4, 1, 12) )
# mat <- model.matrix( ~ a + b, df )
# A <- t(mat) %*% mat # A is SPB

( b <- as.matrix( rnorm(5) ) )

## 1: Compute the Cholesky factorization, \mathbf{A} = \mathbf{LL}^{T} ##
( esky <- chol(A) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
L <- t(esky) ; tL <- esky

round( L %*% tL, 14 ) == round( A, 14 )

## 2: Solve \mathbf{L}\bm{\upsilon} = \mathbf{b} (forward substitution) ##
diag(L) == diag(tL)
diag(L) != 0

( upsi <- numeric( nrow(b) ) )

upsi[1] <- b[1] / L[1, 1]


for (i in 2:length(b))
  upsi[i] <- ( 1 / L[i, i] ) *
  ( b[i] - sum( L[i, 1:(i-1)] * upsi[1:(i-1)] ) ) ; upsi

## 3: Solve \mathbf{L}^{T}\mathbf{x} = \bm{\upsilon} ================== ##
##    Back substitution =============================================== ##
( x <- numeric( length(upsi) ) )

x[5] <- upsi[5] / L[5, 5]

for (i in (length(upsi)-1):1)
  x[i] <-
  ( 1 / L[i, i] ) *
  ( upsi[i] - sum( L[(i+1):length(upsi), i] * x[(i+1):length(upsi)] ) )
## 4: Return \mathbf{x} ============================================== ##
x

cbind( A %*% x, b )
round( A %*% x, 14 ) == round( b, 14 )

# Algorithm 2.2 =========================================================
##  Solving \mathbf{AX} = \mathbf{B} where \mathbf{A} > 0 ============ ##
A

( B <- mvtnorm::rmvnorm(n = 5, mean = rep(0, 7), sigma = diag(7)) )

## 1: Compute the Cholesky factorization, \mathbf{A} = \mathbf{LL}^{T} ##
L
tL

round( L %*% tL, 14 ) == round( A, 14 )

## 2: for j = 1 to k do ============================================== ##
ncol(B) # k

## 3: - Solve \mathbf{L}\bm{\upsilon} = \mathbf{B}_{j} =============== ##
##      (forward substitution) ======================================= ##
## 4: - Solve \mathbf{L}^{T}\mathbf{X}_{j} = \bm{\upsilon} =========== ##
##      Back substitution ============================================ ##
( X <- matrix(0, nrow = nrow(B), ncol = ncol(B)) )

for (j in 1:ncol(B)) {
  
  upsi <- numeric( nrow(B) )
  upsi[1] <- B[1, j] / L[1, 1]
  
  for (i in 2:5)
    upsi[i] <-
    ( 1 / L[i, i] ) * ( B[i, j] - sum( L[i, 1:(i-1)] * upsi[1:(i-1)] ) )
  
  X[5, j] <- upsi[5] / L[5, 5]
  for (i in (nrow(B)-1):1)
    X[i, j] <-
    ( 1 / L[i, i] ) *
    ( upsi[i] - sum( L[(i+1):nrow(B), i] * X[(i+1):nrow(B), j] ) )
}
## 5: end for ======================================================== ##
## 6: Return \mathbf{X} ============================================== ##
X

A %*% X
B
round( A %*% X, 13 ) == round( B, 13 )

## 2.3.2 Unconditional simulation of a GMRF =============================

## Algorithm 2.3 ========================================================
## Sampling \mathbf{x} \sim N(\bm{\mu}, \bm{\Sigma}) ================= ##
( mu <- as.matrix(rnorm(5), nrow = 5) )
( sig <- cov( mvtnorm::rmvnorm(n = 10
                               , mean = rep(0, nrow(mu))
                               , sigma = diag(nrow(mu))) ) ) # sig is SPB

## 1: Compute the Cholesky factorization, ============================ ##
##    \bm{\Sigma} = \tilde{\mathbf{L}}\tilde{\mathbf{L}}^{T} ========= ##
( esky <- chol(sig) ) # this matrix is upper triangular,
                      # but \tilde{L} have to be lower, so,
                      # esky = \tilde{L}^{T}
t(esky) # this matrix is lower triangular,
        # but \tilde{L}^{T} have to be upper, so esky^{T} = \tilde{L}
L <- t(esky) ; tL <- esky

round( L %*% tL, 14 ) == round( sig, 14 )

## 2: Sample \mathbf{z} \sim N(\mathbf{0}, \mathbf{I}) =============== ##
( z <- mvtnorm::rmvnorm(n = nrow(mu)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )

## 3: Compute \bm{\upsilon} = \tilde{\mathbf{L}}\mathbf{z} =========== ##
( upsi <- L %*% z )

## 4: Compute \mathbf{x} = \bm{\mu} + \bm{\upsilon} ================== ##
x <- mu + upsi

## 5: Return \mathbf{x} ============================================== ##
x

## Algorithm 2.4 ======================================================== 
## Sampling \mathbf{x} \sim N(\bm{\mu}, \mathbf{Q}^{-1}) ============= ##
mu
( q <- solve(sig) )

## 1: Compute the Cholesky factorization, ============================ ##
##    \mathbf{Q} = \mathbf{L}\mathbf{L}^{T} ========================== ##
( esky <- chol(q) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
L <- t(esky) ; tL <- esky

round( L %*% tL, 14 ) == round( q, 14 )

## 2: Sample \mathbf{z} \sim N(\mathbf{0}, \mathbf{I}) =============== ##
( z <- mvtnorm::rmvnorm(n = nrow(mu)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )

## 3: Solve \mathbf{L}^{T}\bm{\upsilon} = \mathbf{z} ================= ##
##    Back substitution ============================================== ##
( upsi <- numeric( nrow(mu) ) )

upsi[5] <- z[5] / L[5, 5]

for (i in (nrow(z)-1):1)
  upsi[i] <-
  ( 1 / L[i, i] ) *
  ( z[i] - sum( L[(i+1):nrow(z), i] * upsi[(i+1):nrow(z)] ) ) ; upsi

## 4: Compute \mathbf{x} = \bm{\mu} + \bm{\upsilon} ================== ##
x <- mu + upsi

## 5: Return \mathbf{x} ============================================== ##
x

### ================================================================= ###