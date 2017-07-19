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

all.equal( L %*% tL, A )

## 2: Solve \mathbf{L}\bm{\upsilon} = \mathbf{b} (forward substitution) ##
all.equal( diag(L), diag(tL) )
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
all.equal(A %*% x, b)

# Algorithm 2.2 =========================================================
##  Solving \mathbf{AX} = \mathbf{B} where \mathbf{A} > 0 ============ ##
A

( B <- mvtnorm::rmvnorm(n = 5, mean = rep(0, 7), sigma = diag(7)) )

## 1: Compute the Cholesky factorization, \mathbf{A} = \mathbf{LL}^{T} ##
L
tL

all.equal( L %*% tL, A )

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
all.equal( A %*% X, B )

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

all.equal( L %*% tL, sig )

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

all.equal( L %*% tL, Q )

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

all.equal( L %*% tL, Q )

## 2: Solve \mathbf{L}\bm{\omega} = \mathbf{b} ======================= ##
##    Forward substitution =========================================== ##
# all.equal( diag(L), diag(tL) )
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

all.equal( L %*% tL, Q )

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

all.equal( Lw %*% tLw, W )

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

all.equal( L %*% tL, Q )

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

all.equal( Lw %*% tLw, W )

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

all.equal( L_se %*% tL_se, sig_eps )

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

## 2.4 Numerical methods for sparce matrices ============================
## 2.4.1 Factorizing a sparce matrix ====================================

# Algorithm 2.8 =========================================================
## Cholesky factorization of \mathbf{Q} > 0 ========================== ##
Q

## 1: for j = 1 to n do ============================================== ##
## 2: - \upsilon_{j:n} = Q_{j:n, j} ================================== ##
## 3: - for k = 1 to j - 1 do ======================================== ##
##      \upsilon_{j:n} = \upsilon_{j:n} - L_{j:n, k} L_{j, k} ======== ##
## 4: - L_{j:n, j} = \upsilon_{j:n} / \sqrt{\upsilon_{j}} ============ ##
n <- ncol(Q)

L <- matrix(0, nrow = n, ncol = n)

L[ , 1] <- Q[ , 1] / sqrt( Q[1, 1] )

for (j in 2:n) {
  upsi <- Q[j:n, j]
  
  for (k in 1:(j-1)) upsi <- upsi - L[j:n, k] * L[j, k]
  
  L[j:n, j] <- upsi / sqrt(upsi[1])
}
## 5: end for ======================================================== ##
## 6: Return \mathbf{L} ============================================== ##
L

all.equal( L, t(chol(Q)) )

# Algorithm 2.9 =========================================================
## Band-Cholesky factorization of \mathbf{Q} with bandwidth p ======== ##
# n <- 3 ; p <- 1
# Q <- Matrix::sparseMatrix(i = c(1:2, 1:3, 2:3)
#                           , j = c(1, 1, 2, 2, 2, 3, 3)
#                           , x = c(2, -1, -1, 2, -1, -1, 2)
#                           , dims = c(n, n) ) ; Q

n <- 7 ; p <- 2
Q <- Matrix::sparseMatrix(
  i = c( 1:(p+1), 1:(p+2)
         , rep(5:n, each = 5) + c(-4:0), (n-p-1):n, (n-p):n )
  , j = c( rep(1, p+1), rep(2, p+2)
           , rep((p+1):(n-p), each = 5), rep(n-1, p+2), rep(n, p+1) )
  , x = c( 5, -1, -1, -1, 5, -1, -1
           , rep(c(-1, -1, 5, -1, -1), n-2*p), -1, -1, 5, -1, -1, -1, 5 )
  , dims = c(n, n) ) ; Q

## 1: for j = 1 to n do ============================================== ##
## 2: - \lambda = min{j + p, n} ====================================== ##
## 3: - \upsilon_{j:k} = Q_{j:\lambda, j} ============================ ##
## 4: - for k = max{1, j - p} to j - 1 do ============================ ##
## 5:   - i = min{k + p, n} ========================================== ##
## 6:   - \upsilon_{j:i} = \upsilon_{j:i} - L_{j:i, k} L_{j, k} ====== ##
## 7: - end for ====================================================== ##
## 8: - L_{j:\lambda, j} = \upsilon_{j:\lambda} / \sqrt{\upsilon_{j}}  ##
upsi <- numeric(n)

L <- matrix(0, nrow = n, ncol = n)

for (j in 1:n) {
  l <- min(j + p, n)
  
  upsi[j:l] <- Q[j:l, j]
  
  if (j > 1)
    for (k in max(1, j-p):(j-1)) {
      
      i <- min(k + p, n)
      
      upsi[j:i] <- upsi[j:i] - L[j:i, k] * L[j, k] }
  
  L[j:l, j] <- upsi[j:l] / sqrt(upsi[j])
}
## 9: end for ======================================================== ##
## 10: Return \mathbf{L} ============================================= ##
L

all.equal( L, t(chol(Q)), check.attributes = FALSE )

## 2.6 Stationary GMRFs =================================================
## 2.6.3 GMRFs with circulant precision matrices ========================

# Algorithm 2.10 ========================================================
## Sampling a zero mean GMRF with block-circulant precision ========== ##

cm <- function(x) { # cm: circulant matrix
  n <- length(x)
  suppressWarnings(
    matrix(x[matrix(1:n, n+1, n+1, byrow = TRUE)[c(1, n:2), 1:n]], n, n)
  )
} ; cm(letters[1:4])
cols <- function(ms, n) { # ms: matrices (list of)
  if (n == 0) return(ms)
  c( tail(ms, n), head(ms, -n) )
}
rcols <- function(n, ms) do.call(rbind, cols(ms, n))
bcm <- function(ms) { # bcm: Block-Circulant Matrix
  n <- length(ms)
  do.call( cbind, lapply(0:(n-1), rcols, ms) )
}
### Three different circulating matrix of dimension 3 x 3
( ms <- list( cm(c(5, 2, 2)), cm(c(2, 0, 0)), cm(1:3) ) )
( th <- bcm(ms) )
chol(th) # Showing that the resulting matrix is SPB

## 1: Sample \mathbf{z}, where Re(z_{ij}) \overset{iid}{\sim} N(0, 1) and
##    Im(z_{ij}) \overset{iid}{\sim} N(0, 1) ========================= ##

( z <- rnorm(nrow(th)) + rnorm(nrow(th)) * 1i )

## 2: Compute the (real) eigenvalues, ================================ ##
##    \bm{\Lambda} = \sqrt(nN) DFT2(\theta) ========================== ##
## 3: \bm{\upsilon} = DFT2( ( \bm{\Lambda} \textcircled{e}(-frac{1}{2}) )
##                          \odot \mathbf{z} ) ======================= ##
## 4: \mathbf{x} = Re(\bm{\upsilon})
## 5: Return mathbf{x} =============================================== ##

### ================================================================= ###