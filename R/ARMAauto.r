ARMAauto <- function(ar = NULL, ma = NULL, lag.max)
{

	##########################################################################
	#
	#	ARMAauto
	# 	    Copyright (C) 2017  Tucker McElroy
	#
	#    This program is free software: you can redistribute it and/or modify
	#    it under the terms of the GNU General Public License as published by
	#    the Free Software Foundation, either version 3 of the License, or
	#    (at your option) any later version.
	#
	#    This program is distributed in the hope that it will be useful,
	#    but WITHOUT ANY WARRANTY; without even the implied warranty of
	#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#    GNU General Public License for more details.
	#
	#    You should have received a copy of the GNU General Public License
	#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	#
	############################################################################

	################# Documentation #####################################
	#
	#	Purpose: compute the autocovariance function of an ARMA process
	#	Background: function computes autocovariances of ARMA (p,q) from lag zero
	#		to lag.max, with inputs ar and ma.  Format:
	#		(1 - ar[1]z ... - ar[p]z^p) X_t = (1 + ma[1]z ...+ ma[q]z^q) WN
	#           For absent AR or MA portions, pass in NULL
	#	Inputs:
	#		ar: numeric vector of AR coefficients
	#		ma: numeric vector of MA coefficients
	#		lag.max: Largest autocovariance lag required
	#	Outputs:
	#		autocovariances at lags 0 through lag.max
	#	Requires: polymult
	#
	##########################################################################################

	p <- length(ar)
	q <- length(ma)
	gamMA <- polymult(c(1,ma),rev(c(1,ma)))
	gamMA <- gamMA[(q+1):(2*q+1)]

	if (p > 0)
	{
		Amat <- matrix(0,nrow=(p+1),ncol=(2*p+1))
		for(i in 1:(p+1))
		{
			Amat[i,i:(i+p)] <- c(-1*rev(ar),1)
		}
		Amat <- cbind(Amat[,(p+1)],as.matrix(Amat[,(p+2):(2*p+1)]) +
				t(matrix(apply(t(matrix(Amat[,1:p],p+1,p)),2,rev),p,p+1)))
		Bmat <- matrix(0,nrow=(q+1),ncol=(p+q+1))
		for(i in 1:(q+1))
		{
			Bmat[i,i:(i+p)] <- c(-1*rev(ar),1)
		}
		Bmat <- t(matrix(apply(t(Bmat),2,rev),p+q+1,q+1))
		Bmat <- matrix(apply(Bmat,2,rev),q+1,p+q+1)
		Bmat <- Bmat[,1:(q+1)]
		Binv <- solve(Bmat)

		gamMix <- Binv %*% gamMA
		if (p <= q) gamMix <- matrix(gamMix[1:(p+1),],p+1,1) else
		{ gamMix <- matrix(c(gamMix,rep(0,(p-q))),p+1,1) }
		gamARMA <- solve(Amat) %*% gamMix
	} else gamARMA <- gamMA[1]

	gamMA <- as.vector(gamMA)
	if (lag.max <= q) gamMA <- gamMA[1:(lag.max+1)] else gamMA <- c(gamMA,rep(0,(lag.max-q)))
	gamARMA <- as.vector(gamARMA)
	if (lag.max <= p) gamARMA <- gamARMA[1:(lag.max+1)] else {
		for(k in 1:(lag.max-p))
		{
			len <- length(gamARMA)
			acf <- gamMA[p+1+k]
			if (p > 0) acf <- acf + sum(ar*rev(gamARMA[(len-p+1):len]))
			gamARMA <- c(gamARMA,acf)
		}
	}
	return(gamARMA)
}
