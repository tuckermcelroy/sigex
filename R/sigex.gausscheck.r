sigex.gausscheck <- function(resids)
{

	##########################################################################
	#
	#	sigex.gausscheck
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
	#	Purpose: wrapper for Shapiro-Wilks test of normality
	#	Background:
	#		model-fitting is an entropy maximizing transformation of
	#		the data, producing residuals that should resemble
	#		Gaussian white noise.  We can test whether the 
	#		marginal distribution is normal with Shapiro-Wilks.
	#	Inputs:
	#		resids: a T x N matrix of residuals
	#	Outputs:
	#		tests: output of shapiro.tests routine of R
	#
	####################################################################
	
	x <- t(resids)
	N <- dim(x)[1]
	T <- dim(x)[2]
#	if(T > 5000) { T <- 5000 }

	tests <- NULL
	for(i in 1:N)
	{
		tests <- c(tests,shapiro.test(resids[1:T,i])$p.value)
	}

	return(tests)
}


