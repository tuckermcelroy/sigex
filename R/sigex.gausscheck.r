sigex.gausscheck <- function(resid)
{

	#################################
	#	sigex.gausscheck
	#		by Tucker McElroy
	#
	#	wrapper for Shapiro-Wilks test of normality
	#
	#########################################

	x <- t(resid)
	N <- dim(x)[1]
	T <- dim(x)[2]
	if(T > 5000) { T <- 5000 }

	tests <- NULL
	for(i in 1:N)
	{
		tests <- c(tests,shapiro.test(resid[1:T,i])$p.value)
	}

	return(tests)
}


