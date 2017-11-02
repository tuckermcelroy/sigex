sigex.param2gcd <- function(L.psi,N,vrank)
{

	###########################
	#	sigex.param2gcd
	#		by Tucker McElroy
	#
	#   Takes parameter vector psiL and returns unit lower
	#	triangular matrix of dimension N x length(vrank) with
	#	entries given by L.psi
	#
	##############################


	L.mat <- 1
	if(N > 1) {
	fill <- NULL
	ind <- 0
	for(j in vrank)
	{
		pad <- NULL
		if(j < N) pad <- L.psi[(ind+1):(ind + N-j)]
		new <- c(rep(0,j-1),1,pad)
		fill <- c(fill,new)
		ind <- ind + N - j
	}
	L.mat <- matrix(fill,nrow=N,ncol=length(vrank)) 
	}

	return(L.mat)
}

