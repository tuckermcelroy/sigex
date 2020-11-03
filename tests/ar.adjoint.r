ar.adjoint <- function(poly.array)
{

	#########################################
	#	ar.adjoint by Tucker McElroy
	#
	#	given matrix polynomial array poly.array
	#	 of dimension N,N,p+1 (usual plus convention),
	#	 where the leading coefficient matrix need not
	#	 be identity, returns adjoint matrix polynomial
	#	 as array of dimension N,N,p(N-1)+1
	#	This version assumes p >= 1
	#	Returns list with adjoint and determinant polynomial
	#########################################

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


	
	

