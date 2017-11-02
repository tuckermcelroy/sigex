sigex.glr <- function(data,psi.nested,psi.nesting,mdl.nested,mdl.nesting)
{

	########################################
	#
	#	SigEx GLR
	#		Tucker McElroy
	#
	#	Computes the difference of -2*loglik for
	#		two models, the nested lik minus nesting lik
	#		(this can be applied to non-nested models,
	#		but the distribution won't be chi^2
	#
	#########################################

	glr <- sigex.lik(psi.nested,mdl.nested,data) - 
		sigex.lik(psi.nesting,mdl.nesting,data)
	dof <- length(psi.nesting) - length(psi.nested)

	return(c(glr,dof))
}
