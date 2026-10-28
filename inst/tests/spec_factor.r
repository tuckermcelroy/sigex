polymult <- function(a,b)
{
  bb <- c(b,rep(0,length(a)-1))
  B <- toeplitz(bb)
  B[lower.tri(B)] <- 0
  aa <- rev(c(a,rep(0,length(b)-1)))
  prod <- B %*% matrix(aa,length(aa),1)
  return(rev(prod[,1]))
}

spec_factor <- function(poly,eps=10^{-3})
{
  ############################################################
  # spec_factor by Tucker McElroy, October 2026
  # Purpose: compute spectral factorization of a given
  #   symmetric 2-sided polynomial
  # Inputs:
  #   poly: input vector of coefficients, of length 2*p.2+1
  #     NOTE: must be symmetric (poly == rev(poly)) and 
  #           must be non-negative definite ( FT(poly) >= 0 )
  #   eps: minimum phase distance between root clusters
  # Outputs:
  #   factor: output vector of coefficients, of length p.2+1,
  #     such that poly == polymult(factor,rev(factor))
  # Requires: polymult
  
  p <- length(poly)-1  # needs p even
  theta <- 1
  prod <- poly[p+1]
  roots <- polyroot(poly)
#  roots <- round(roots,digits=8)
#  roots <- pracma::roots(poly)

  # need to make sure that all roots have -pi < Arg <= pi, 
  # by adding 2 pi to any argument equal to -pi
  roots.arg <- Arg(roots)
  ind.mpi <- which(abs(roots.arg + pi) < eps)
  if(length(ind.mpi)>0) { roots.arg[ind.mpi] <- roots.arg[ind.mpi] + 2*pi }
  
  # sort roots by Arg
  roots.ind <- order(roots.arg)
  roots <- roots[roots.ind]
  roots.mat <- matrix(roots,nrow=2)
  p.2 <- dim(roots.mat)[2] # should be p/2
  
  # find root clusters and build polynomial factor
#  cluster.begs <- c(0,which(diff(Arg(roots))>eps)/2) +1
  cluster.begs <- c(0,which(diff(sort(roots.arg))>eps)/2) +1
  for(j in 1:length(cluster.begs))
  {
    if(j == length(cluster.begs))
    {
      cluster.end <- p.2
    } else
    {
      cluster.end <- cluster.begs[j+1]-1
    }
    # Get all roots of the cluster
    cluster.roots <- roots.mat[,cluster.begs[j]:cluster.end]
    # sort cluster roots by modulus, should have an even number of roots
    cluster.roots <- cluster.roots[order(Mod(cluster.roots))] 
    # associate first half of cluster roots with the factor
    for(k in 1:(length(cluster.roots)/2))
    {
      theta <- polymult(theta,c(1,-cluster.roots[k]))
      prod <- prod/(-cluster.roots[k])
    }
  }

  factor <- Re(theta)*sqrt(Re(prod))
  return(factor)
}

## Test cases

library(pracma)

eps <- 10^{-3}

### (1+B^2)(1-B)(1+rho^2 B^2)
rho <- .3
factor <- polymult(polymult(c(1,0,1),c(1,-1)),c(1,0,rho^2))
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))
 
### (1 - B^2)
factor <- c(1,0,-1)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

### (1 + B)^3
factor <- c(1,3,3,1)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

### 1 + B^3
factor <- c(1,0,0,1)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

### 1-B^4
factor <- c(1,0,0,0,-1)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

### (1 + B)(1 - B)^2
factor <- c(1,-1,-1,1)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

### (1+B)^6
factor <- choose(6, 0:6)
poly <- polymult(factor,rev(factor))
print(rbind(factor,spec_factor(poly)))

#####################################################
## SCRAP

# Separate out unit roots from stable/explosive roots
inds.unit <- which(abs(Mod(roots)-1)<eps)
if(length(inds.unit)>0)
{
  if(length(inds.unit)<p)
  {
    roots.stable <- roots[-inds.unit]
    roots.unit <- roots[inds.unit]
  }
  else
  {
    roots.stable <- NULL
    roots.unit <- roots   
  }
} else
{
  roots.stable <- roots
  roots.unit <- NULL
}

# Get contribution from stable/explosive roots
if(length(inds.unit)<p)
{
  p.stable <- p - length(inds.unit)
  for(i in 1:p.stable)
  {
    if (Mod(roots.stable[i]) < 1)
    {
      theta <- polymult(theta,c(1,-roots[i]))
      prod <- prod/(-roots[i])
    } 
  }
}
 
