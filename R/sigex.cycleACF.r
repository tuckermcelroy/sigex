sigex.cycleACF <- function(cycle.order,rho,omega,delta,maxlag)
{

	########################################
	#
	#	sigex.cycleACF
	#		by Tucker McElroy
	#
	#	gets acf of a Balanced cycle that has been differenced
	#		with polynomial delta, with persistency rho
	#		and frequency omega.
	#	Gives acf for lags 0 through maxlag-1
	#
	######################

	d.delta <- length(delta)
	cycle.acf <- c(1,rep(0,maxlag+d.delta-1))

	# JL - I got an ERROR on line 20 so I commented it out
	#choose(cycle.order-1

	out <- filter(rev(filter(testacf,c(1,-1),method="convolution",sides=2)),c(1,-1),method="convolution",sides=2)


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


