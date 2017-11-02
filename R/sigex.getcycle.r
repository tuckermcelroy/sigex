sigex.getcycle <- function(cycle.order,rho,omega)
{

	########################################
	#
	#	sigex.getcycle
	#		by Tucker McElroy	
	#
	#	get cycle AR and MA polynomials (minus and plus conventions)
	#		corresponding to the Butterworth cycle
	#
	######################

	cycle.AR <- 1
	cycle.MA <- 1
	if(cycle.order > 0)
	{
		if(omega==0) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymul(cycle.AR,c(1,-rho))
		} }
		if(omega==1) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymul(cycle.AR,c(1,rho))
		} }
		if((omega > 0) & (omega < 1)) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymul(cycle.AR,c(1,-2*rho*cos(pi*omega),rho^2))
			cycle.MA <- polymul(cycle.MA,c(1,-1*rho*cos(pi*omega)))
		} 
		#cycle.acf <- polymul(cycle.MA,rev(cycle.MA))
		}
	}

	return(list(cycle.AR,cycle.MA))
}


