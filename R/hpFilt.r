hpFilt <- function(q,n)
{
	# Written by Tucker McElroy of U.S. Census Bureau, July 2006
	# Input Signal-Noise ratio q
	# Total number of coefficients is 2n + 1
	# absZ is the absolute value of zeta, the AR polynomial root
	 absZ <- (sqrt(q) + sqrt(q+16) + sqrt(2*q + 2*sqrt(q)*sqrt(q+16)))/4
	# c is a normalizing constant
	 c <- q/(absZ^2)
	# theta is the angular portion of zeta 
	 theta <- atan( sqrt(2*q + 2*sqrt(q)*sqrt(q+16))/4)
	# phi's are the HP filter coefficients 
	 phi <- matrix(0,n+1,1)
	 for(j in 0:n)
	 {
	   phi[j+1] <- 2*c*absZ^(4-j)*sin(theta)*(absZ^2*sin((j+1)*theta)-sin((j-1)*theta))/
			((absZ^2 - 1)*(1-2*absZ^2*cos(2*theta)+absZ^4)*(1-cos(2*theta)))
	 }
	 filter_coef<-phi
	# Other output that may be of interest
	#  phi1 is the first AR coefficient
	 phi1 <- 2*(sqrt(q) - sqrt(q+16))/(sqrt(q) + sqrt(q+16) + sqrt(2*q + 2*sqrt(q)*sqrt(q+16)))
	# phi2 is the second AR coefficient
	 phi2 <- (sqrt(q) + sqrt(q+16) - sqrt(2*q + 2*sqrt(q)*sqrt(q+16)))/
		 (sqrt(q) + sqrt(q+16) + sqrt(2*q + 2*sqrt(q)*sqrt(q+16)))
	 ma_model<-c(c,phi1,phi2)
	# For the roots, note that absZ*cos(theta) + i*absZ*sin(theta) gives the one in
	#  the first quadrant, and absZ*cos(theta) - i*absZ*sin(theta) is in the fourth quadrant
#	 print(c)
#	 print(phi1)
#	 print(phi2)
#	 print(absZ)
#	 print(theta)
#	 print(absZ*cos(theta))
#	 print(absZ*sin(theta))
 	 return(list(filter_coef=filter_coef,ma_model=ma_model))	
}