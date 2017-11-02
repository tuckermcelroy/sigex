sigex.delta <- function(mdl,omits)
{
	#################################
	#   sigex.delta
	#	by Tucker McElroy	
	#
	#	Multiplies all the delta polynomials for each series, 
	#		except the ones corresponding to indices of "omits";
	#		has no effect if that component is already 
	#			omitted from a series
	#
	#################################

	prod <- 1
	for(i in 1:length(mdl[[3]]))
	{
		polyn <- mdl[[3]][[i]]
		if (i %in% omits) polyn <- 1
		prod <- polymul(prod,polyn)		
	}
	return(prod)
}



