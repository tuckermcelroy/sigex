sigex.coordesc <- function(psi.init,fix.lik,psi.fix,mdl,data.ts,flag,lower.bound,upper.bound)
{

	##########################################################################
	#
	#	sigex.coordesc
	# 	    Copyright (C) 2019  Tucker McElroy
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
	#	Purpose: optimize via coordinate descent
	#	Background:	
	#		Wrapper for a sigex likelihood, intended as an
	#		alternative to BFGS or CG methods	
	#	Inputs:
	#		fix.lik: likelihood with multiple arguments, whose
	#			first argument should be for the parameter psi
	#		psi.init: initial value of the parameter vector
	#	Outputs:
	#	Requires: 
	####################################################################
  
	p <- length(psi.init)
	tol <- 10^(-10)
	do.it <- TRUE
	psi.old <- psi.init
	psi.new <- psi.old
	my.coord <- 1
	while(do.it)
	{	  	
		onedim.lik <- function(psi1,psi.rest,mdl,data.ts,my.coord)
		{
			psi.est <- rep(0,p)
			psi.est[-my.coord] <- psi.rest
			psi.est[my.coord] <- psi1	
			out <- fix.lik(psi.est,psi.fix,mdl,data.ts,flag)
			return(out)
		}
		psi.rest <- psi.old[-my.coord]
		onedim.fit <- optimize(onedim.lik,c(lower.bound[my.coord],upper.bound[my.coord]),
			psi.rest=psi.rest,mdl=mdl,data.ts=data.ts,flag=flag)
		psi.new[my.coord] <- onedim.fit$par
		my.coord <- my.coord + 1
		if(my.coord == p+1) 
		{
			my.coord <- 1
			psi.old <- psi.new
		}
		do.it <- sum((psi.new - psi.old)^2) > tol			
	}

	return(psi.new)
}





