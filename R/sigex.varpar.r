sigex.varpar <- function(psi,var.order,N,delta,debug=FALSE)
{

	##########################################################################
	#
	#	sigex.varpar.r  
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
	#	Purpose: generates a stable dimension N VAR(p) process
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		an array object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#		Algorithm is that of Roy, McElroy, Linton.
	#	Notes: only for use with N > 1 and p > 0
	#	Inputs:
	#		psi: vector of real numbers, of length p*N^2
	#		var.order: the VAR order p 
	#		N: dimension of the process
	#		delta: p-vector of plus/minus one
	#		debug: a Boolean flag; if true, outputs stability check
	#	Outputs:
	#		param: N x N x var.order array
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

	psi.mat <- matrix(psi,nrow=N*N,ncol=var.order)
	psi.mat <- cbind(psi.mat,rep(0,N*N))
	v.mat <- array(0, dim = c(N, N, var.order))
	Q.mat <- v.mat

	for(j in 1:var.order)
	{
		l.mat <- diag(N)
		l.mat[lower.tri(l.mat)] <- psi.mat[1:choose(N, 2), j]
		d.mat <- diag(exp(psi.mat[(choose(N, 2) + 1):choose(N + 1, 2), j]))
  		v.mat[,,j] <- l.mat %*% d.mat %*% t(l.mat)
		s.mat <- diag(0, N)
		s.mat[lower.tri(s.mat)] <- psi.mat[(choose(N + 1, 2) + 1):(N^2), j]
   		s.mat <- s.mat - t(s.mat)
   		Q.mat[,,j] <- diag(c(delta[j], rep(1, (N - 1)))) %*% 
			(diag(N) - s.mat) %*% solve(diag(N) + s.mat)
  	}

	l.mat <- diag(N)
	l.mat[lower.tri(l.mat)] <- psi.mat[1:choose(N, 2),(var.order + 1)]
	d.mat <- diag(exp(psi.mat[(choose(N,2)+1):choose(N + 1, 2),(var.order + 1)]))

  	param <- array(0, dim = c(N,N,var.order))
	u <- array(0, dim = c(N,N,(var.order + 1)))
  	u[,,1] <- diag(N) + apply(v.mat, c(1, 2), sum)
  	dd <- u[,,1]
	u[,,2] <- sqrtm(v.mat[,,1]) %*% Q.mat[,,1] %*% sqrtm(dd)
	tkappa <- t(u[,,2])
	txi <- u[,,2]
  	bigu <- u[,,1]
	if (var.order > 1) {
		for(j in 3:(var.order + 1)){
          		dd <- u[,,1] - tkappa %*% solve(bigu) %*% t(tkappa)
	      	u[,,j] <- (txi %*% solve(bigu) %*% t(tkappa) + 
      	            sqrtm(v.mat[,,(j-1)]) %*% Q.mat[,,(j-1)] %*% sqrtm(dd) )
      		bigu <- rbind(cbind(bigu,t(tkappa)), cbind(tkappa,u[,,1]))
	  	      tkappa <- cbind(t(u[,,j]), tkappa)
	   		txi <- cbind(txi, u[,,j])
    		}
  	}
  	phi.mat <- txi %*% solve(bigu)
	for(j in 1:var.order) { param[,,j] <- phi.mat[,(N * (j - 1) + 1):(N * j)] }
   
	# stability check
	if(debug) {
		phis.var <- matrix(param[,,1:var.order],nrow=N)
		if(var.order==1) { trans.mat <- phis.var } else {
			trans.mat <- rbind(phis.var,diag(N*var.order)[1:(N*(var.order-1)),]) }
		if(max(Mod(eigen(trans.mat)$value)) < 1) { print("Stable!") }
	}

	return(param)

}
