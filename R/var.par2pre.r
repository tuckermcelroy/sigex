var.par2pre <- function(param)
{
  
  ##########################################################################
  #
  #	var.par2pre.r  
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
  #	Purpose: generates a stable dimension N VAR(p) process
  #	Background:	
  #		param is the name for the model parameters entered into 
  #		an array object with a more intuitive structure, whereas
  #		psi refers to a vector of real numbers containing all
  #		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
  #		Algorithm is that of Ansley and Kohn (1986)
  #	Notes: only for use with N > 1 and p > 0
  #	Inputs:
  #		param: N x N x var.order array,
  #		  where var.order is the VAR order p and N is the dimension of the process
  #	Outputs:
  #		psi: vector of real numbers, of length p*N^2  
  # Requires: VARMAauto.r
  #
  ####################################################################
  
  sqrtm <- function(A) { return(t(chol(A))) }

  N <- dim(param)[1]
  var.order <- dim(param)[3]
  acvf <- VARMAauto(param,NULL,diag(N),var.order)
  
  Sigma <- acvf[,,1]
  S.mat <- sqrtm(Sigma)
  Sstar.mat <- S.mat
  sigma.array <- array(0,c(N,N,var.order+1))
  sigmastar.array <- array(0,c(N,N,var.order+1))
  sigma.array[,,1] <- Sigma
  sigmastar.array[,,1] <- Sigma
  phi.array <- list()
  phistar.array <- list()
  pacf.array <- array(0,c(N,N,var.order))
  for(k in 1:var.order)
  {
    phi.array[[k]] <- array(0,c(N,N,k))
    phistar.array[[k]] <- array(0,c(N,N,k))
    new.phi <- acvf[,,k+1] %*% solve(sigmastar.array[,,k])
    new.phistar <- t(acvf[,,k+1]) %*% solve(sigma.array[,,k])
    if(k > 1)
    {
      for(j in 1:(k-1))
      {
        new.phi <- new.phi - phi.array[[k-1]][,,j] %*% acvf[,,k-j+1] %*% solve(sigmastar.array[,,k])
        new.phistar <- new.phistar - phistar.array[[k-1]][,,j] %*% t(acvf[,,k-j+1]) %*% solve(sigma.array[,,k])
      }
    }
    phi.array[[k]][,,k] <- new.phi
    phistar.array[[k]][,,k] <- new.phistar
    if(k > 1)
    {
      for(j in 1:(k-1))
      {
        phi.array[[k]][,,j] <- phi.array[[k-1]][,,j] - phi.array[[k]][,,k] %*% phistar.array[[k-1]][,,k-j]
        phistar.array[[k]][,,j] <- phistar.array[[k-1]][,,j] - phistar.array[[k]][,,k] %*% phi.array[[k-1]][,,k-j]
      }
    }
    pacf.array[,,k] <- solve(S.mat) %*% phi.array[[k]][,,k] %*% Sstar.mat
 
    new.sigma <- acvf[,,1]
    new.sigmastar <- acvf[,,1]
    for(j in 1:k)
    {
        new.sigma <- new.sigma - phi.array[[k]][,,j] %*% t(acvf[,,j+1])
        new.sigmastar <- new.sigmastar - phistar.array[[k]][,,j] %*% acvf[,,j+1]
    }
    sigma.array[,,k+1] <- new.sigma
    sigmastar.array[,,k+1] <- new.sigmastar

    if(k < var.order)
    {
      S.mat <- sqrtm(sigma.array[,,k+1])
      Sstar.mat <- sqrtm(sigmastar.array[,,k+1])
    }
  }

  psi.array <- array(0,c(N,N,var.order))
  for(k in 1:var.order)
  { 
    psi.array[,,k] <- solve(sqrtm(diag(N) - pacf.array[,,k] %*% 
                              t(pacf.array[,,k]))) %*% pacf.array[,,k]
  }
  psi.mat <- matrix(psi.array,c(N,N*var.order))
  psi <- matrix(psi.mat,ncol=1)    
  
  return(psi)
}
