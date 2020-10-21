var2.par2pre <- function(param)
{
  
  ##########################################################################
  #
  #	var2.par2pre.r  
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
  #		Algorithm is that of Roy, McElroy, and Linton (2019)
  #	Notes: only for use with N > 1 and p > 0
  #	Inputs:
  #		param: N x N x p array,
  #		  where p is the VAR order and N is the dimension of the process
  #	Outputs:
  #		psi: vector of real numbers, of length p*N^2  
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
  
  N <- dim(param)[1]
  p <- dim(param)[3]
  u <- VARMAauto(param,NULL,diag(N),p)
  A <- array(0, dim = c(N,N,p))
  dd <- u[,,1]
  ddinv <- solve(dd)
  A[,,1] <- u[,,2] %*% msqrt(dd)$invsqrt
  T1 <- 0*A[,,1]
  txi <- tkappa <- T2 <- ddinv <- NULL
  if (p > 1){
    for(j in 2:p){
      if(j == 2){uinv <- solve(u[,,1]) }
      else{uinv <- rbind(cbind((uinv + T2%*%t(T1)), -T2), cbind(-t(T2),ddinv))  }
      tkappa <- cbind(t(u[,,j]),tkappa)
      txi <- cbind(txi, u[,,j])
      T1 <- uinv %*% t(tkappa)
      dd <- u[,,1] - tkappa %*% T1
      ddinv <- solve(dd)
      T2 <- T1 %*% ddinv
      A[,,j] <- (u[,,(j+1)] - txi %*% uinv %*% t(tkappa)) %*% msqrt(dd)$invsqrt
    }
  }
  pre <- c(A)
  return(pre)
}

