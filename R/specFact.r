specFact <- function(poly)
{
	p <- length(poly)-1
	roots <- polyroot(poly)
	theta <- 1
	prod <- poly[p+1]
	toggle <- 1
	for(i in 1:p)
	{
		if (Mod(roots[i]) < 1)
		{
			theta <- polymult(theta,c(1,-roots[i]))
			prod <- prod/(-roots[i])
		} else {
		if (Mod(roots[i]) <= 1) 
		{
			if(Arg(roots[i]) < 0)
			{
				theta <- polymult(theta,c(1,-roots[i]))
				prod <- prod/(-roots[i])
			} else {
			if((Arg(roots[i]) == 0) && (toggle == 1))
			{	
				theta <- polymult(theta,c(1,-roots[i]))
				prod <- prod/(-roots[i])
				toggle <- -1*toggle
			} }
		} }
	}
	theta <- Re(theta)*sqrt(Re(prod))
	return(theta)
}
		  