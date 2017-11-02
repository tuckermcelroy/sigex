sigex.meaninit <- function(mdl,data,d)
{
	#################################
	#   sigex.meaninit
	#	by Tucker McElroy	
	#
	#	Adds trend regressors to an existing model
	#		to all of the   series; always use this
	#		function, before calling sigex.reg  
	#		if no trend component exists, then a polynomial
	#		 trend of order d is utilized; otherwise,
	#		 d is obtained from the trend differencing order
	#	Model is described via mdl, a list object:
	#		mdl[[1]] is mdlK gives ranks of components, from vrank
	#		mdl[[2]] is mdlType gives model types for components
	#		mdl[[3]] is mdlDiff (a list)
	#			gives delta differencing polynomials, from delta
	#		mdl[[4]] is list of regressors by series
	#		mdl[[5]] is list of bounds for rho, omega
	#
	#	x is NxT   time series
	#
	#################################

mdlK <- mdl[[1]]
mdlType <- mdl[[2]]
mdlDiff <- mdl[[3]]
mdlReg <- mdl[[4]]
mdlBounds <- mdl[[5]]
 
if(length(sigex.whichtrend(mdl))==1) {
delta <- mdl[[3]][[sigex.whichtrend(mdl)]]
d <- length(delta) - 1 - sum(abs(Arg(polyroot(delta))) > 10^(-8))
#d <- length(delta) - 1 
}
x <- t(data)
T <- dim(x)[2]
N <- dim(x)[1]
delta <- sigex.delta(mdl,0)
for(series in 1:N)
{ 
	mdlReg[[series]] <- ts(as.matrix(seq(1,T)^d),start=start(data),
		frequency=frequency(data),names="Trend") 
	d.inc <- d
	while(d.inc > 0)
	{
		d.inc <- d.inc-1
		reg <- ts(as.matrix(seq(1,T)^d.inc),start=start(data),
			frequency=frequency(data),names="Trend")
		reg.diff <- filter(reg,delta,method="convolution",
			sides=1)[length(delta):T] 
		if(sum(reg.diff^2) > 10^(-8)) 
		{
			mdlReg[[series]] <- ts(cbind(mdlReg[[series]],reg),
				start=start(reg),frequency=frequency(reg),
				names=c(colnames(mdlReg[[series]]),colnames(reg)))
		}
	}
}
mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg,bounds = mdlBounds)
return(mdl)
}

