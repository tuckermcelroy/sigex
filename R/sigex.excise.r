sigex.excise <- function(data,param,mdl,excise.comp,excise.series)
{

	#################################
	#	sigex.excise
	#		by Tucker McElroy
	#
	#	Constructs a parameter-constrained model
	#		where certain components for particular
	#		series are excised.
	#	excise.comp indicates which components are excised
	#	excise.series indicates for which series it is excised
	#		these must have the same length
	#
	#	Output is new mdl, and corresponding par (with fixed values
	#		flagged with 1i)
	#
	#########################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	mdl.new <- mdl
	par.new <- param
	for(j in unique(excise.comp))
	{
		excise.j <- excise.series[excise.comp==j]
		ranks <- mdl[[1]][[j]]
		new.ranks <- setdiff(ranks,excise.j)
		mdl.new[[1]][[j]] <- new.ranks
	
		par.new[[1]][[j]] <- par.new[[1]][[j]][,new.ranks]
		if(length(new.ranks)==1) par.new[[1]][[j]] <- matrix(par.new[[1]][[j]],ncol=1)
		for(k in excise.j) { par.new[[1]][[j]][k,] <- 0+1i }
		par.new[[2]][[j]] <- par.new[[2]][[j]][new.ranks] 
	}

	#	psi.new <- sigex.input(par.new,mdl.new)

	return(list(mdl.new,par.new))
}
	




	