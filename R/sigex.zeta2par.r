sigex.zeta2par <- function(zeta,mdlType,N,bounds)
{

	##########################################################################
	#
	#	sigex.zeta2par
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
	#	Purpose: transform zeta to param
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#	Notes: this is a functional inverse to sigex.par2zeta
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Inputs:
	#		zeta: see background.  This is the portion of psi
	#			corresponding t.s. models, such as cycles
	#		mdlType: this is a component of mdl (the specified sigex model),
	#			 cited as mdl[[2]]
	#		N: dimension of the time series
	#		bounds: gives bounds for rho and omega, cycle parameters in zeta
	#			rho lies in (bounds[1],bounds[2])
	#			omega lies in (bounds[3],bounds[4])
	#	Outputs:
	#		zeta.par: see background.  This is a portion of the full 
	#			param list, corresponding to param[[3]]
	#	Requires: sigex.param2gcd
	#
	####################################################################

	# defaults
	low.rho <- bounds[1]
	upp.rho <- bounds[2]
	low.omega <- bounds[3]
	upp.omega <- bounds[4]
		
	if(mdlType %in% c("wn","canonWN")) { zeta.par <- NULL }
	if(mdlType == "AR1") 
	{ 
		phi1 <- low.rho + (upp.rho - low.rho)*exp(zeta[1])/(1+exp(zeta[1]))
		zeta.par <- phi1
	} 
	if(mdlType %in% c("MA1","canonMA1")) 
	{ 
		theta1 <- (exp(zeta[1])-1)/(exp(zeta[1])+1)
		zeta.par <- theta1
	} 
	if(mdlType %in% c("cycleBW1","cycleBW2","cycleBW3","cycleBW4","cycleBW5",
		"cycleBW6","cycleBW7","cycleBW8","cycleBW9","cycleBW10",
		"canonCycleBW1","canonCycleBW2","canonCycleBW3","canonCycleBW4",
		"canonCycleBW5","canonCycleBW6","canonCycleBW7","canonCycleBW8",
		"canonCycleBW9","canonCycleBW10","cycleBAL1","cycleBAL2","cycleBAL3",
		"cycleBAL4","cycleBAL5","cycleBAL6","cycleBAL7","cycleBAL8",
		"cycleBAL9","cycleBAL10","canonCycleBAL1","canonCycleBAL2",
		"canonCycleBAL3","canonCycleBAL4","canonCycleBAL5","canonCycleBAL6",
		"canonCycleBAL7","canonCycleBAL8","canonCycleBAL9","canonCycleBAL10"))
	{
		rho <- low.rho + (upp.rho - low.rho)*exp(zeta[1])/(1+exp(zeta[1]))
		omega <- low.omega + (upp.omega - low.omega)*exp(zeta[2])/(1+exp(zeta[2]))
		zeta.par <- c(rho,omega)		
	}
#	if(mdlType == "VAR1") 
#	{ 
#		phi1.mat <- psi[(ind+1):(ind+N^2)]
#		subparam <- matrix(phi1.mat,ncol=N)	
#	}
	if(mdlType == "ARMA22")
	{
		phis <- zeta[1:2]
		phis[1] <- (exp(phis[1])-1)/(exp(phis[1])+1)
		phis[2] <- (exp(phis[2])+2*phis[1]^2-1)/(exp(phis[2])+1)			
		phis <- c(phis[1]*(1-phis[2]),(phis[2]-phis[1]^2))/(1-phis[1]^2)
		thetas <- zeta[3:4]
		thetas[1] <- (exp(thetas[1])-1)/(exp(thetas[1])+1)
		thetas[2] <- (exp(thetas[2])+2*thetas[1]^2-1)/(exp(thetas[2])+1)
		thetas <- c(thetas[1]*(1-thetas[2]),(thetas[2]-thetas[1]^2))/(1-thetas[1]^2)
		zeta.par <- c(phis,-1*thetas)
	}
	if(mdlType == "damped")
	{
		L.dim <- as.integer(N*(N-1)/2)
		L.psi <- zeta[1:L.dim]
		D.psi <- zeta[(L.dim+1):(L.dim+N)]
		L.mat <- sigex.param2gcd(L.psi,N,seq(1,N))
		phi <- zeta[L.dim+N+1]
		zeta.par <- list(L.mat,D.psi,phi)
	} 

	return(zeta.par)
}
