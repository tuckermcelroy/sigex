sigex.renderpd <- function(L.mat,D.mat,thresh)
{
	
	########################################
	#	SigEx renderpd
	#		by Tucker McElroy
	#
	#	Takes given covariance matrix in L D L' format,
	#		and returns a modified D such that
	#		condition numbers of new cov matrix with
 	#		the same L are bounded below by thresh.
	#	Note: if run on a reduced rank matrix, pad out
	#		D.mat with -Inf and corresponding columns
	#		of L.mat should be unit vector.
	#	Dimension >= 2, because first condition number is zero.
	#
	##################################	

	N <- dim(L.mat)[1]
	D.new <- matrix(exp(D.mat[1]),1,1)
	for(i in 2:N)
	{
		val <- (matrix(L.mat[i,1:(i-1)],nrow=1) %*% D.new %*% matrix(L.mat[i,1:(i-1)],ncol=1))/(exp(-thresh) - 1)
		d.new <- max(exp(D.mat[i]),val)
		D.new <- diag(c(diag(D.new),d.new))
	}
	return(log(diag(D.new)))
}
