m <- 30
p <- 2
phi <- array(0, dim = c(m, m, p))
for(i in 1:p){
  phi[,,i] <- matrix(rnorm(m^2, 0, .002), m, m)
}

bigphi <- rbind(matrix(phi, m, m*p), cbind(diag(m*(p-1)), matrix(0, m*(p-1), m)))
max(abs(eigen(bigphi)$values))
sigma = diag(m)

g <- VARMAauto(phi, NULL, sigma, p)
