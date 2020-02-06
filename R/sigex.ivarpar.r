sigex.ivarpar <- function(param)
{

	##########################################################################
	#
	#	sigex.ivarpar.r  
	# 	    Copyright (C) 2019  Tucker McElroy
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
	#	Purpose: computes a stable dimension N VAR(p) process
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		an array object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#		Algorithm is that of Roy, McElroy, Linton.
	#	Notes: only for use with N > 1 and p > 0
	#	Inputs:
	#		param: N x N x p array
	#	Outputs:
	#		psi: vector of real numbers, of length p*N^2
	#		delta: p-vector of plus/minus one
	#	Requires: getGCD
	#
	####################################################################

sqrtm <- function(A) {
	# Function returns the square root of the matrix
	# Input: 
	#   A: m x m non-negative definite matrix
	# Output:
	#   m x m matrix
	return(eigen(A)$vectors %*% diag(sqrt(eigen(A)$values)) %*% t(eigen(A)$vectors))
}

	N <- dim(param)[1]
  p <- dim(param)[3]
	phi.mat <- matrix(param,N,N*p)
	phi.comp <- rbind(phi.mat,cbind(diag(N*(p-1)),matrix(0,N*(p-1),N)))
	v.mat <- array(0,dim=c(N,N,p))
	q.mat <- v.mat
	u.mat = array(0,dim=c(N,N,(p+1)))
	c.mat <- u.mat
	d.mat <- u.mat
	if(p==1)
	{
		U.bmat <- matrix(solve((diag(N^2*p^2) - phi.comp %x% phi.comp),
			as.vector(diag(N)),tol=1e-40),N*p,N*p)
		v.mat[,,1] = U.bmat - diag(N)
            q.mat[,,1] = solve(sqrtm(v.mat[,,1])) %*% phi.mat %*% sqrtm(U.bmat)
            u.mat[,,1] = U.bmat
            u.mat[,,2] = phi.mat %*% U.bmat
	} else 
	{
		U.bmat <- matrix(solve((diag(N^2*p^2) - phi.comp %x% phi.comp),
			as.vector(diag(c(rep(1,N),rep(0,N*(p-1))))),tol=1e-40),N*p,N*p)
	      for(j in 1:p) { u.mat[,,j] <- U.bmat[1:N,((j-1)*N+1):(j*N)] }
		u.mat[,,(p+1)] <- phi.mat %*% U.bmat[1:(N*p),(N*(p-1)+1):(N*p)]
		d.mat[,,1] <- u.mat[,,1]
    		c.mat[,,1] <- u.mat[,,1]
		bigu <- u.mat[,,1]
		tkappa <- t(u.mat[,,2])
		txi <- u.mat[,,2]
		tt <- u.mat[,,2]
		c.mat[,,2] <- u.mat[,,1] - txi %*% solve(bigu) %*% t(txi)
		d.mat[,,2] <- u.mat[,,1] - tkappa %*% solve(bigu) %*% t(tkappa)
		v.mat[,,1] <- c.mat[,,1] - c.mat[,,2]
		q.mat[,,1] <- solve(sqrtm(v.mat[,,1])) %*% u.mat[,,2] %*% solve(sqrtm(d.mat[,,1]))
		for(j in 2:p) 
		{
			tt <- u.mat[,,(j+1)] - txi %*% solve(bigu) %*% t(tkappa)
                  bigu <- rbind(cbind(bigu,t(tkappa)),cbind(tkappa,u.mat[,,1]))
                  tkappa <- cbind(t(u.mat[,,(j+1)]),tkappa)
                  txi <- cbind(txi,u.mat[,,(j+1)]) 
                  d.mat[,,(j+1)] <- u.mat[,,1] - tkappa %*% solve(bigu) %*% t(tkappa) 
                  c.mat[,,(j+1)] <- u.mat[,,1] - txi %*% solve(bigu) %*% t(txi)
                  v.mat[,,j] <- c.mat[,,j] - c.mat[,,(j+1)]
                  q.mat[,,j] <- solve(sqrtm(v.mat[,,j])) %*% tt %*% solve(sqrtm(d.mat[,,j])) 
		}
	}

	psi <- matrix(0,nrow=N^2,ncol=p)
	delta <- rep(0,p)
	for(j in 1:p)
	{
		out <- getGCD(v.mat[,,j],N)
		psi[1:choose(N+1,2),j] <- c(out[[1]][lower.tri(out[[1]])],log(out[[2]]))
		delta[j] <- sign(det(q.mat[,,j]))
		s.mat <- 2 * solve(diag(N) + diag(c(delta[j],rep(1,(N-1)))) %*% q.mat[, , j]) - diag(N)
		psi[(choose(N+1,2)+1):N^2,j] <- s.mat[lower.tri(s.mat)]
	}
	psi <- as.vector(matrix(psi,ncol=1))

	return(list(psi,delta))
}


