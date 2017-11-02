sigex.reduce <- function(data,param,flag,mdl,thresh,modelflag)
{

	#################################
	#	sigex.reduce
	#		by Tucker McElroy
	#
	#	Computes a reduced rank model, by
	#		replacing small Schur comps with omission when modelflag = TRUE
	#		if modelflag = FALSE, we keep model the same and replace small
	#			Schur comps by method involving exp(thresh)
	#	Output is new mdl, and corresponding par 
	#
	#########################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	
	psi.red <- sigex.par2psi(param,flag,mdl)
	log.conds <- log(sigex.conditions(data,psi.red,mdl))

	if(modelflag) {
	par.red <- param
	mdl.red <- NULL
	for(j in 1:length(mdl[[3]])) {
		ranks <- seq(1,N)[log.conds[j,] > thresh]
		mdl.red <- sigex.add(mdl.red,ranks,mdl[[2]][j],mdl[[3]][[j]],mdl[[5]][[j]])
		par.red[[1]][[j]] <- as.matrix(param[[1]][[j]][,ranks])
		par.red[[2]][[j]] <- param[[2]][[j]][ranks]
	}
	mdl.red <- sigex.meaninit(mdl.red,data,0)
	mdl.red[[4]] <- mdl[[4]] } else {
	par.red <- param
	mdl.red <- mdl
	for(j in 1:length(mdl[[3]])) 
	{
		L.mat <- par.red[[1]][[j]]
		D.psi <- par.red[[2]][[j]]
		D.new <- sigex.renderpd(L.mat,D.psi,thresh)
		par.red[[2]][[j]] <- D.new
#		cov.mat <- L.mat %*% diag(exp(D.psi),nrow=length(D.psi)) %*% t(L.mat)
#		ranks <- seq(1,N)[log.conds[j,] < thresh]
#		if(length(ranks) > 0) {
#		par.red[[2]][[j]][ranks] <- rep(thresh,length(ranks)) + log(diag(cov.mat))[ranks] }
	} }

	return(list(mdl.red,par.red))
}

