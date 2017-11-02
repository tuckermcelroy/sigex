sigex.mvar2uvar <- function(data,param,mdl)
{


	###############################
	#   sigex.mvar2uvar
	#	by Tucker McElroy	
	#
	#	Transforms multivariate parameters to implied univariate form
	#		param must be in format yielded by sigex.default
	#     Output is same format as param
	#
	#################################	

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	mdl.uni <- mdl
	for(i in 1:length(mdl.uni[[1]])) { mdl.uni[[1]][[i]] <- seq(1,N) }
	default.param <- sigex.default(mdl.uni,x)
	univ.param <- param
	for(i in 1:length(mdl[[3]]))
	{
		L.psi <- param[[1]][[i]]
		D.psi <- param[[2]][[i]]
		k.order <- length(D.psi)
	 	S.mat <- L.psi %*% diag(exp(D.psi),nrow=k.order) %*% t(L.psi)
		univ.param[[2]][[i]] <- log(diag(S.mat))				
		univ.param[[1]][[i]] <- default.param[[1]][[i]]
	}

	return(list(mdl.uni,univ.param))
}
