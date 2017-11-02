sigex.add <- function(mdl,vrank,type,delta,bounds)
{
	#################################
	#   sigex.add
	#	by Tucker McElroy	
	#
	#	Adds new latent components to an existing model
	#		to any of the component series of index in "series"
	#	Model is described via mdl, a list object:
	#		mdl[[1]] is mdlK gives ranks of components, from vrank
	#		mdl[[2]] is mdlType gives model types for components
	#		mdl[[3]] is mdlDiff (a list)
	#			gives delta differencing polynomials, from delta
	#		mdl[[4]] is list of regressors by series
	#		mdl[[5]] is list of bounds for rho, omega
	#
	#################################

mdlK <- mdl[[1]]
mdlType <- mdl[[2]]
mdlDiff <- mdl[[3]]
mdlReg <- mdl[[4]]
mdlBounds <- mdl[[5]]

rank.null <- NULL
if(length(vrank)==1) rank.null <- 0
mdlK[[length(mdlK)+1]] <- c(rank.null,vrank)
if(length(vrank)==1) mdlK[[length(mdlK)]] <- mdlK[[length(mdlK)]][-1]
mdlType <- c(mdlType,type)
delta.null <- NULL
if(length(delta)==1) delta.null <- 0
mdlDiff[[length(mdlDiff)+1]] <- c(delta.null,delta)
if(length(delta)==1) mdlDiff[[length(mdlDiff)]] <- mdlDiff[[length(mdlDiff)]][-1]
mdlBounds[[length(mdlBounds)+1]] <- bounds 

mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg,bounds = mdlBounds)
return(mdl)
}

