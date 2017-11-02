sigex.blocktoep <- function(x.array)
{
	
	#################################
	#	sigex.blocktoep
	#
	#	takes x.array of dimension N,N,H
	#		and generates block Toeplitz lower triangular 
	#		array of dimension N,H,N,H
	#
	##############################

N <- dim(x.array)[1]
H <- dim(x.array)[3]

x.toep <- array(0,c(N,H,N,H))
for(j in 1:H)
{
	for(k in 1:H)
	{
		if(j >= k) { x.toep[,j,,k] <- x.array[,,j-k+1] }
	}
}

return(x.toep)
}
