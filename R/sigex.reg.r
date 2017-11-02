sigex.reg <- function(mdl,series,reg)
{
	#################################
	#   sigex.reg
	#	by Tucker McElroy	
	#
	#	Adds new regressors to an existing model
	#		to any of the   series  
	#	Model is described via mdl, a list object:
	#		mdl[[1]] is mdlK gives ranks of components, from vrank
	#		mdl[[2]] is mdlType gives model types for components
	#		mdl[[3]] is mdlDiff (a list)
	#			gives delta differencing polynomials, from delta
	#		mdl[[4]] is list of regressors by series
	#	series is between 1 and N, the series for which regressors (reg)
	#		are being added.  
	#	reg is a vector of time series regressors, of length T
	#		(presumes that correct start and end dates for series 
	#		 of length T are known, and regressors are constructed accordingly)
	#	Regressors are entered; if magnitude of differences is small, they are deemed
	#		to be in the null space of the diff ops, and hence are omitted
	#		to avoid misidentification.  
	#
	#################################

mdlK <- mdl[[1]]
mdlType <- mdl[[2]]
mdlDiff <- mdl[[3]]
mdlReg <- mdl[[4]]

T <- length(reg)
delta <-  sigex.delta(mdl,0)
reg.diff <- filter(reg,delta,method="convolution",sides=1)[length(delta):T] 
if(sum(reg.diff^2) > 10^(-8)) 
{
	mdlReg[[series]] <- ts(cbind(mdlReg[[series]],reg),start=start(reg),
		frequency=frequency(reg),names=c(colnames(mdlReg[[series]]),colnames(reg)))
}
mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg)

return(mdl)
}

