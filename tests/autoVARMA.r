autoVARMA <- function(phi,theta,sigma,grid=1000,maxlag)
{
  
  ##########################################################################
  #
  #	autoVARMA
  # 	    Copyright (C) 2020  Tucker McElroy
  #
  #    This program is free software: you can redistribute it and/or modify
  #    it under the terms of the GNU General Public License as published by
  #    the Free Software Foundation, either version 3 of the License, or
  #    (at your option) any later version.
  #
  #    This program is distributed in the hope that it will be useful,
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #    GNU General Public License for more details.
  #
  #    You should have received a copy of the GNU General Public License
  #    along with this program.  If not, see <https://www.gnu.org/licenses/>.
  #
  ############################################################################
  
  ################# Documentation #####################################
  #
  #	Purpose: computes autocovariances of VARMA usng frequency domain
  #	Background: function computes autocovariances of VARMA (p,q) from lag zero
  #		to maxlag, with array inputs phi and theta.  VARMA equation:
  #	(1 - phi[1]B ... - phi[p]B^p) X_t = (1 + theta[1]B ...+ theta[q]B^q) WN_t
  #	Note: for absent VAR or VMA portions, pass in NULL
  #	Inputs:
  #		phi: array of dimension m x m x p of VAR coefficients, e.g.,
  #			phi <- array(cbind(phi1,phi2,...,phip),c(m,m,p))
  #		theta: array of dimension m x m x q of VMA coefficients, e.g.,
  #			theta <- array(cbind(theta1,theta2,...,thetaq),c(m,m,q))
  #		sigma: m x m covariance matrix of white noise
  #   grid: Riemann mesh size
  #	Outputs:
  #		autocovariances at lags 0 through maxlag, as array of dimension m x m x (maxlag+1)
  #
  ####################################################################

  polymult <- function(a,b) 
  {
    bb <- c(b,rep(0,length(a)-1))
    B <- toeplitz(bb)
    B[lower.tri(B)] <- 0
    aa <- rev(c(a,rep(0,length(b)-1)))
    prod <- B %*% matrix(aa,length(aa),1)
    return(rev(prod[,1]))
  }
  
  polymulMat <- function(amat,bmat)
  {
    p <- dim(amat)[3]-1
    q <- dim(bmat)[3]-1
    N <- dim(amat)[2]
    
    r <- p+q
    bmat.pad <- array(0,c(N,N,r+1))
    for(i in 1:(q+1)) { bmat.pad[,,i] <- bmat[,,i] }
    cmat <- array(0,c(N,N,r+1))
    cmat[,,1] <- amat[,,1] %*% bmat.pad[,,1]
    for(j in 2:(r+1))
    {
      cmat[,,j] <- amat[,,1] %*% bmat.pad[,,j]
      for(k in 1:min(p,j-1))
      { cmat[,,j] <- cmat[,,j] + amat[,,k+1] %*% bmat.pad[,,j-k] }
    }
    
    return(cmat)
  }

  ar.adjoint <- function(poly.array)
  {
    p <- dim(poly.array)[3]-1
    N <- dim(poly.array)[1]
    
    poly.0 <- poly.array[,,1]
    poly.mat <- matrix(poly.array[,,2:(p+1)],c(N,p*N))
    poly.coefs <- -1*solve(poly.0) %*% poly.mat
    
    if(p==1) { companion.mat <- poly.coefs } else {
      companion.mat <- diag(p*N)[1:((p-1)*N),]	
      companion.mat <- rbind(poly.coefs,companion.mat) }
    poly.evals <- eigen(companion.mat)$values
    ar.poly <- det(poly.0)
    for(j in 1:(p*N))
    { ar.poly <- polymult(ar.poly,c(1,-poly.evals[j])) }
    ar.poly <- Re(ar.poly)
    
    r <- p*(N-1)+1
    adj.array <- array(0,c(N,N,r)) 
    adj.array[,,1] <- ar.poly[1]*solve(poly.0)
    for(j in 2:r)
    {
      adj.array[,,j] <- ar.poly[j]*solve(poly.0)
      for(k in 1:min(p,j-1))
      {
        adj.array[,,j] <- adj.array[,,j] - solve(poly.0) %*% 
          poly.array[,,k+1] %*% adj.array[,,j-k]
      }
    }
    
    return(list(adj.array,ar.poly))
  }
  
  N <- dim(sigma)[1]
  p <- dim(phi)[3]
  q <- dim(theta)[3]
  phi.long <- array(cbind(diag(N),matrix(phi,c(N,N*p))),c(N,N,p+1))
  out <- ar.adjoint(phi.long)
  phi.adjoint <- out[[1]]
  phi.det <- out[[2]]
  
  freqs <- seq(-1*grid,grid)*pi/grid
  len <- length(freqs)
  lambdas <- exp(-1i*freqs)
  phi.adjointz <- array(0,c(N,N*len))
  for(k in 1:dim(phi.adjoint)[3])
  {
    phi.adjointz <- phi.adjointz + t(lambdas^(k-1)) %x% phi.adjoint[,,k]
  }
  phi.adjointz <- array(phi.adjointz,c(N,N,len))
  phi.detz <- rep(0,len)
  for(k in 1:length(phi.det))
  {
    phi.detz <- phi.detz + lambdas^(k-1) * phi.det[k]
  }
  
  gamma <- array(0,c(N,N,maxlag+1))
  for(h in 0:maxlag)
  {
    val <- do.call(cbind,lapply(seq(1,len),function(i) phi.adjointz[,,i] %*% 
                sigma %*% t(Conj(phi.adjointz[,,i])) * Mod(phi.detz[i])^{-2} * lambdas[i]^{-h}))
    val <- len^{-1}*val %*% (rep(1,len) %x% diag(N))
    gamma[,,h+1] <- Re(val)
  }  
  
  return(gamma)
  
}