sigex.portmanteau <- function(resids,lag,nump)
{

	##########################################################################
	#
	#	sigex.portmanteau
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
	#	Purpose: computes the portmanteau statistic for residuals
	#	Background:
	#		model-fitting is an entropy maximizing transformation of
	#		the data, producing residuals that should resemble
	#		Gaussian white noise.  We can test whether there is  
	#		serial autocorrelation via the portmanteau statistic
	#	Inputs:
	#		resids: a T x N matrix of residuals
	#		lag: number of autocorrelation lags used in the statistic
	#		nump: number of parameters that are estimated
	#	Outputs:
	#		port: value of the portmanteau test statistic
	#		pval: p-value of port, based on chi square with nump 
	#			degrees of freedom 
	#
	####################################################################
 
	x <- t(resids)
	N <- dim(x)[1]
	T <- dim(x)[2]

	dof <- lag*N^2 - nump
	if(dof <= 0) { lag <- nump + 1; dof <- 1 }
	acf.sample <- acf(t(x),type="covariance",lag.max=lag,plot=FALSE)$acf
	varinv <- solve(acf.sample[1,,])
	port <- 0
	for(h in 1:lag)
	{
		port <- port + sum(diag(acf.sample[h+1,,] %*% varinv %*% 
			t(acf.sample[h+1,,]) %*% varinv))
	}
	port <- T*port
	pval <- 1-pchisq(port,df= dof)
	
	return(c(port,pval))
}



