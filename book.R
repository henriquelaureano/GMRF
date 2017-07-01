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
  
  for (i in 2:nrow(B))
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
( Q <- solve(sig) )

## 1: Compute the Cholesky factorization, ============================ ##
##    \mathbf{Q} = \mathbf{L}\mathbf{L}^{T} ========================== ##
( esky <- chol(Q) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
L <- t(esky) ; tL <- esky

round( L %*% tL, 14 ) == round( Q, 14 )

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

## Algorithm 2.5 ========================================================
## Sampling \mathbf{x} \sim N_{C}(\mathbf{b}, \mathbf{Q}) ============ ## 
b
Q

## 1: Compute the Cholesky factorization, \mathbf{Q} = \mathbf{LL}^{T} ##
L # t(chol(q))
tL

round( L %*% tL, 13 ) == round( Q, 13 )

## 2: Solve \mathbf{L}\bm{\omega} = \mathbf{b} ======================= ##
##    Forward substitution =========================================== ##
# diag(L) == diag(tL)
# diag(L) != 0

( om <- numeric( nrow(b) ) )

om[1] <- b[1] / L[1, 1]

for (i in 2:length(b))
  om[i] <- ( 1 / L[i, i] ) *
  ( b[i] - sum( L[i, 1:(i-1)] * om[1:(i-1)] ) ) ; om

## 3: Solve \mathbf{L}^{T}\bm{\mu} = \bm{\omega} ===================== ##
##    Back substitution ============================================== ##
( mu <- numeric( length(om) ) )

mu[5] <- om[5] / L[5, 5]

for (i in (length(om)-1):1)
  mu[i] <-
  ( 1 / L[i, i] ) *
  ( om[i] - sum( L[(i+1):length(om), i] * mu[(i+1):length(om)] ) ) ; mu

## 4: Sample \mathbf{z} \sim N(\mathbf{0}, \mathbf{I}) =============== ##
( z <- mvtnorm::rmvnorm(n = length(mu)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )

## 5: Solve \mathbf{L}^{T}\bm{\upsilon} = \mathbf{z} ================= ##
( upsi <- numeric( nrow(z) ) )

upsi[5] <- z[5] / L[5, 5]

for (i in (length(upsi)-1):1)
  upsi[i] <-
  ( 1 / L[i, i] ) *
  ( z[i] - sum( L[(i+1):nrow(z), i] * upsi[(i+1):nrow(z)] ) ) ; upsi

## 6: Compute \mathbf{x} = \bm{\mu} + \bm{\upsilon} ================== ##
x <- mu + upsi

## 7: Return \mathbf{x} ============================================== ##
x

## 2.3.3 Conditional simulation of a GMRF ===============================

## Algorithm 2.6 ========================================================
## Sampling \mathbf{x}|\mathbf{Ax} = \mathbf{e} ====================== ##
## where \mathbf{x} \sim N(\bm{\mu}, \mathbf{Q}^{-1}) ================ ##
mu
Q
( A <- mvtnorm::rmvnorm(n = 4, mean = rep(0, 5), sigma = diag(5)) )
( e <- rnorm( nrow(A) ) )

## 1: Compute the Cholesky factorization, \mathbf{Q} = \mathbf{LL}^{T} ##
L # t(chol(Q))
tL

round( L %*% tL, 14 ) == round( Q, 14 )

## 2: Sample \mathbf{z} \sim N(\mathbf{0}, \mathbf{I}) =============== ##
( z <- mvtnorm::rmvnorm(n = length(mu)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )

## 3: Solve \mathbf{L}^{T}\bm{\upsilon} = \mathbf{z} ================= ##
##    Back substitution ============================================== ##
( upsi <- numeric( nrow(z) ) )

upsi[5] <- z[5] / L[5, 5]

for (i in (nrow(z)-1):1)
  upsi[i] <-
  ( 1 / L[i, i] ) *
  ( z[i] - sum( L[(i+1):nrow(z), i] * upsi[(i+1):nrow(z)] ) ) ; upsi

## 4: Compute \mathbf{x} = \bm{\mu} + \bm{\upsilon} ================== ##
( x <- mu + upsi)

## 5: Compute \mathbf{V}_{n \times k} = \mathbf{Q}^{-1}\mathbf{A}^{T}  ##
##    using Algorithm 2.2 using \mathbf{L} from step 1 =============== ##

## ========= \mathbf{V} = \mathbf{Q}^{-1}\mathbf{A}^{T} \Rightarrow == ##
## \mathbf{Q}\mathbf{V} = \mathbf{A}^{T} ============================= ##
( V <- matrix( 0, nrow = nrow(Q), ncol = ncol(t(A)) ) )

 for (j in 1:ncol(t(A))) {

  upsi <- numeric( nrow(t(A)) )
  upsi[1] <- t(A)[1, j] / L[1, 1]

  for (i in 2:nrow(t(A)))
    upsi[i] <-
    ( 1 / L[i, i] ) *
    ( t(A)[i, j] - sum( L[i, 1:(i-1)] * upsi[1:(i-1)] ) )

  V[5, j] <- upsi[5] / L[5, 5]
  for (i in (nrow(t(A))-1):1)
    V[i, j] <-
    ( 1 / L[i, i] ) *
    ( upsi[i] - sum( L[(i+1):nrow(t(A)), i] * V[(i+1):nrow(t(A)), j] ) )
} ; V

## 6: Compute \mathbf{W}_{k \times k} = \mathbf{AV} ================== ##
( W <- A %*% V )

## 7: Compute \mathbf{U}_{k \times n} = \mathbf{W}^{-1}\mathbf{V}^{T}  ##
##    using Algorithm 2.2 ============================================ ##

## ========= \mathbf{U} = \mathbf{W}^{-1}\mathbf{V}^{T} \Rightarrow == ##
## \mathbf{W}\mathbf{U} = \mathbf{V}^{T} ============================= ##

( esky <- chol(W) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
Lw <- t(esky) ; tLw <- esky

round( Lw %*% tLw, 14 ) == round( W, 14 )

( U <- matrix( 0, nrow = nrow(W), ncol = ncol(t(V)) ) )

for (j in 1:ncol(t(V))) {
  
  upsi <- numeric( nrow(t(V)) )
  upsi[1] <- t(V)[1, j] / Lw[1, 1]
  
  for (i in 2:nrow(t(V)))
    upsi[i] <-
    ( 1 / Lw[i, i] ) *
    ( t(V)[i, j] - sum( Lw[i, 1:(i-1)] * upsi[1:(i-1)] ) )
  
  U[4, j] <- upsi[4] / Lw[4, 4]
  for (i in (nrow(t(V))-1):1)
    U[i, j] <-
    ( 1 / Lw[i, i] ) *
    ( upsi[i] - sum( Lw[(i+1):nrow(t(V)), i] * U[(i+1):nrow(t(V)), j] ) )
} ; U

## 8: Compute \mathbf{c} = \mathbf{Ax} - \mathbf{e} ================== ##
( c <- A %*% x - e )

## 9: Compute \mathbf{x}^{*} = \mathbf{x} - \mathbf{U}^{T}\mathbf{c} = ##
xs <- x - t(U) %*% c

## 10: Return \mathbf{x}^{*} ========================================= ##
xs

## Algorithm 2.7 ========================================================
## Sampling \pi(\mathbf{x}|\mathbf{e}) where ========================= ##
## \mathbf{x} \sim N(\bm{\mu}, \mathbf{Q}^{-1}) and ================== ##
## \mathbf{e}|\mathbf{x} \sim N(\mathbf{Ax}, \bm{\Sigma_{\epsilon}}) = ##
mu
Q
e
( A <- mvtnorm::rmvnorm(n = 4, mean = rep(0, 5), sigma = diag(5)) )
( sig_eps <- cov(
  mvtnorm::rmvnorm(n = 10
                   , mean = rep(0, nrow(A))
                   , sigma = diag(nrow(A)))
) ) # sig_eps is SPB

## 1: Compute the Cholesky factorization, \mathbf{Q} = \mathbf{LL}^{T} ##
L # t(chol(Q))
tL

round( L %*% tL, 14 ) == round( Q, 14 )

## 2: Sample \mathbf{z} \sim N(\mathbf{0}, \mathbf{I}) =============== ##
( z <- mvtnorm::rmvnorm(n = length(mu)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )

## 3: Solve \mathbf{L}^{T}\bm{\upsilon} = \mathbf{z} ================= ##
##    Back substitution ============================================== ##
( upsi <- numeric( nrow(z) ) )

upsi[5] <- z[5] / L[5, 5]

for (i in (nrow(z)-1):1)
  upsi[i] <-
  ( 1 / L[i, i] ) *
  ( z[i] - sum( L[(i+1):nrow(z), i] * upsi[(i+1):nrow(z)] ) ) ; upsi

## 4: Compute \mathbf{x} = \bm{\mu} + \bm{\upsilon} ================== ##
( x <- mu + upsi)

## 5: Compute \mathbf{V}_{n \times k} = \mathbf{Q}^{-1}\mathbf{A}^{T}  ##
##    using Algorithm 2.2 using \mathbf{L} from step 1 =============== ##

## ========= \mathbf{V} = \mathbf{Q}^{-1}\mathbf{A}^{T} \Rightarrow == ##
## \mathbf{Q}\mathbf{V} = \mathbf{A}^{T} ============================= ##
( V <- matrix( 0, nrow = nrow(Q), ncol = ncol(t(A)) ) )

for (j in 1:ncol(t(A))) {
  
  upsi <- numeric( nrow(t(A)) )
  upsi[1] <- t(A)[1, j] / L[1, 1]
  
  for (i in 2:nrow(t(A)))
    upsi[i] <-
    ( 1 / L[i, i] ) *
    ( t(A)[i, j] - sum( L[i, 1:(i-1)] * upsi[1:(i-1)] ) )
  
  V[5, j] <- upsi[5] / L[5, 5]
  for (i in (nrow(t(A))-1):1)
    V[i, j] <-
    ( 1 / L[i, i] ) *
    ( upsi[i] - sum( L[(i+1):nrow(t(A)), i] * V[(i+1):nrow(t(A)), j] ) )
} ; V

## 6: Compute ======================================================== ##
##    \mathbf{W}_{k \times k} = \mathbf{AV} + \bm{\Sigma_{\epsilon}} = ##
( W <- A %*% V + sig_eps )

## 7: Compute \mathbf{U}_{k \times n} = \mathbf{W}^{-1}\mathbf{V}^{T}  ##
##    using Algorithm 2.2 ============================================ ##

## ========= \mathbf{U} = \mathbf{W}^{-1}\mathbf{V}^{T} \Rightarrow == ##
## \mathbf{W}\mathbf{U} = \mathbf{V}^{T} ============================= ##

( esky <- chol(W) ) # this matrix is upper triangular,
                    # but L have to be lower, so, esky = L^{T}
t(esky) # this matrix is lower triangular,
        # but L^{T} have to be upper, so esky^{T} = L
Lw <- t(esky) ; tLw <- esky

round( Lw %*% tLw, 14 ) == round( W, 14 )

( U <- matrix( 0, nrow = nrow(W), ncol = ncol(t(V)) ) )

for (j in 1:ncol(t(V))) {
  
  upsi <- numeric( nrow(t(V)) )
  upsi[1] <- t(V)[1, j] / Lw[1, 1]
  
  for (i in 2:nrow(t(V)))
    upsi[i] <-
    ( 1 / Lw[i, i] ) *
    ( t(V)[i, j] - sum( Lw[i, 1:(i-1)] * upsi[1:(i-1)] ) )
  
  U[4, j] <- upsi[4] / Lw[4, 4]
  for (i in (nrow(t(V))-1):1)
    U[i, j] <-
    ( 1 / Lw[i, i] ) *
    ( upsi[i] - sum( Lw[(i+1):nrow(t(V)), i] * U[(i+1):nrow(t(V)), j] ) )
} ; U

## 8: Sample \bm{\epsilon} \sim N(\mathbf{e}, \bm{\Sigma_{epsilon}}) = ##
##    using Algorithm 2.3 ============================================ ##

L_se <- t( chol(sig_eps) ) ; tL_se <- chol(sig_eps)

round( L_se %*% tL_se, 15 ) == round( sig_eps, 15 )

( z <- mvtnorm::rmvnorm(n = length(e)
                        , mean = rep(0, 1)
                        , sigma = diag(1)) )
( upsi <- L_se %*% z )

( eps <- e + upsi )

## 9: Compute \mathbf{c} = \mathbf{Ax} - \bm{\epsilon} =============== ##
( c <- A %*% x - eps )

## 10: Compute \mathbf{x}^{*} = \mathbf{x} - \mathbf{U}^{T}\mathbf{c}  ##
xs <- x - t(U) %*% c

## 11: Return \mathbf{x}^{*} ========================================= ##
xs

### ================================================================= ###