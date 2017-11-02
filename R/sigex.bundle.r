sigex.bundle <- function(data,transform,mdl,psi)
{

	#################################
	#	SigEx bundle
	#		by Tucker McElroy
	#	Creates a list object, that bundles
	#		data, transform, span, model, and 
	#		fitted parameters into an "analysis"
	#
	###########################################

	analysis <- list(data = data,transform = transform,model = mdl,psi = psi)

	return(analysis)
}

