var.pre2par <- function(psi,var.order,N,debug=FALSE)
{
  
  ##########################################################################
  #
  #	var.pre2par.r  
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
  #		psi: vector of real numbers, of length p*N^2  
  #		var.order: the VAR order p 
  #		N: dimension of the process
  #		debug: a Boolean flag; if true, outputs stability check
  #	Outputs:
  #		param: N x N x var.order array
  #
  ####################################################################
  
  sqrtm <- function(A) {return(t(chol(A))) }
  
  psi.mat <- matrix(psi,nrow=N,ncol=N*var.order)
  psi.array <- array(psi.mat,c(N,N,var.order))
  pacf.array <- array(0,c(N,N,var.order))
  for(k in 1:var.order)
  { 
    pacf.array[,,k] <- solve(sqrtm(diag(N) + psi.array[,,k] %*% 
                                t(psi.array[,,k]))) %*% psi.array[,,k] 
  }
  
  sqrt.array <- array(0,c(N,N,var.order+1))
  sqrt.array[,,var.order+1] <- diag(N)
  for(k in var.order:1)
  { 
    S.mat <- sqrtm(sqrt.array[,,k+1]) %*% 
      solve(sqrtm(diag(N)-pacf.array[,,k] %*% t(pacf.array[,,k]))) 
    sqrt.array[,,k] <- S.mat %*% t(S.mat)
  }
  Sigma <- sqrt.array[,,1]

  sigma.array <- array(0,c(N,N,var.order+1))
  sigmastar.array <- array(0,c(N,N,var.order+1))
  sigma.array[,,1] <- Sigma
  sigmastar.array[,,1] <- Sigma
  phi.array <- list()
  phistar.array <- list()
  Sstar.mat <- S.mat
  for(k in 1:var.order)
  {
    phi.array[[k]] <- array(0,c(N,N,k))
    phistar.array[[k]] <- array(0,c(N,N,k))
    phi.array[[k]][,,k] <- S.mat %*% pacf.array[,,k] %*% solve(Sstar.mat)
    phistar.array[[k]][,,k] <- Sstar.mat %*% t(pacf.array[,,k]) %*% solve(S.mat)
    if(k > 1)
    {
      for(j in 1:(k-1))
      {
        phi.array[[k]][,,j] <- phi.array[[k-1]][,,j] - phi.array[[k]][,,k] %*% phistar.array[[k-1]][,,k-j]
        phistar.array[[k]][,,j] <- phistar.array[[k-1]][,,j] - phistar.array[[k]][,,k] %*% phi.array[[k-1]][,,k-j]
      }
    }
    sigma.array[,,k+1] <- sigma.array[,,k] - phi.array[[k]][,,k] %*% sigmastar.array[,,k] %*% t(phi.array[[k]][,,k])
    sigmastar.array[,,k+1] <- sigmastar.array[,,k] - phistar.array[[k]][,,k] %*% sigma.array[,,k] %*% t(phistar.array[[k]][,,k])    
    S.mat <- sqrtm(sigma.array[,,k+1])
    Sstar.mat <- sqrtm(sigmastar.array[,,k+1])
  }
  
  # check  sigma.array[,,var.order+1] = diag(N)
  param <- phi.array[[var.order]]

  # stability check
  if(debug) {
    phis.var <- matrix(param[,,1:var.order],nrow=N)
    if(var.order==1) { trans.mat <- phis.var } else {
      trans.mat <- rbind(phis.var,diag(N*var.order)[1:(N*(var.order-1)),]) }
    if(max(Mod(eigen(trans.mat)$value)) < 1) { print("Stable!") }
  }
  
  return(param)
}
