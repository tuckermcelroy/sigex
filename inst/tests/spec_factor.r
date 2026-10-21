spec_factor <- function(poly,eps=10^{-3})
{
  # Spectral factorization of symmetric 2-sided polynomial
  
  p <- length(poly)-1  # needs p even
  theta <- 1
  prod <- poly[p+1]
  roots <- polyroot(poly)
#  roots <- round(roots,digits=8)
  print(roots)

  # need to make sure that all roots have 0 <= Arg <= pi, 
  # by conjugating those that have negative argument
  roots <- Re(roots) + 1i*abs(Im(roots))

  # sort roots by Arg
  roots.ind <- order(Arg(roots))
  roots <- roots[roots.ind]
  roots.mat <- matrix(roots,nrow=2)
  p.2 <- dim(roots.mat)[2] # should be p/2
  
  # find root clusters and build polynomial factor
  cluster.begs <- c(0,which(diff(Arg(roots))>eps))/2 +1
  for(j in 1:length(cluster.begs))
  {
    if(j == length(cluster.begs))
    {
      cluster.end <- p.2
    } else
    {
      cluster.end <- cluster.begs[j+1]-1
    }
    # Get all roots of the cluster, but the bottom row must be conjugated?
    cluster.roots <- roots.mat[1,cluster.begs[j]:cluster.end]
    cluster.roots <- c(cluster.root,
                       Conj(roots.mat[2,cluster.begs[j]:cluster.end]))
    
    # sort cluster roots by modulus, should have an even number of roots
    cluster.roots <- cluster.roots[order(Mod(cluster.roots))] 
    # associate first half of cluster roots with the factor
    for(k in 1:(length(cluster.roots)/2))
    {
      theta <- polymult(theta,c(1,-cluster.roots[k]))
      prod <- prod/(-cluster.roots[k])
    }
  }

  theta <- Re(theta)*sqrt(Re(prod))
  return(theta)
}

  
poly <- polymult(c(1,0,1),c(-1,2,-1))
rho <- .9
poly <- polymult(poly,c(rho^2,0,1+rho^4,0,rho^2))
  
  
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
 
