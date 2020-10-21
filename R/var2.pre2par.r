var2.pre2par <- function(psi,p,N)
{

  ##########################################################################
  #
  #	var2.pre2par.r
  # 	    Copyright (C) 2020  Tucker McElroy
  #       (Original code by Anindya Roy, adapted by Tucker McElroy)
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
  #		p: the VAR order p
  #		N: dimension of the process
  #	Outputs:
  #		param: N x N x p array
  # Requires: VARMAauto.r
  #   (the msqrt function is borrowed from MTS package, but this need not be loaded)
  #
  ####################################################################

  msqrt <- function (M)
  {
    if (!is.matrix(M))
      M = as.matrix(M)
    n1 = nrow(M)
    if (n1 == 1) {
      Mh = sqrt(M)
      Mhinv = 1/Mh
    }
    if (n1 > 1) {
      M = (M + t(M))/2
      m1 = eigen(M)
      V = m1$vectors
      eiv = sqrt(m1$values)
      L = diag(eiv)
      Linv = diag(1/eiv)
      Mh = V %*% L %*% t(V)
      Mhinv = V %*% Linv %*% t(V)
    }
    msqrt <- list(mtxsqrt = Mh, invsqrt = Mhinv)
  }

  prephi <- psi[1:(p*N*N)]
  A <- array(prephi, dim = c(N,N,p))
  param <- array(0,dim = c(N,N,p))
  u <- array(0,dim=c(N,N,(p+1)))
  u[,,1] <- diag(N) + matrix(A,N,N*p) %*% t(matrix(A,N,N*p))
  dd <- u[,,1]
  uinv <- solve(u[,,1])
  u[,,2] <- A[,,1] %*% msqrt(dd)$mtxsqrt
  txi <- u[,,2]
  tkappa <- t(txi)
  if (p > 1){
    for(j in 3:(p+1)){
      T1 <- uinv %*% t(tkappa)
      dd <- u[,,1] - tkappa %*% T1
      dd2 <- msqrt(dd)$mtxsqrt
      u[,,j] <-  txi %*% T1 + A[,,(j-1)] %*% dd2
      ddinv <- solve(dd)
      T2 <- T1 %*% ddinv
      uinv <- rbind(cbind((uinv + T2 %*% t(T1)), -T2), cbind(-t(T2),ddinv))
      tkappa <- cbind(t(u[,,j]),tkappa)
      txi <- cbind(txi, u[,,j])
    }
  }
  phi1 <- txi %*% uinv
  for(j in 1:p){
    param[,,j] <- phi1[,(N*(j-1)+1):(N*j)]
  }
  return(param)
}
