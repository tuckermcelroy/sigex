sigex.par2zeta <- function(mdlPar,mdlType)
{

	##########################################################################
	#
	#	sigex.par2zeta
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
	#	Purpose: transform param to zeta
	#	Background:
	#		param is the name for the model parameters entered into
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Notes: this is a functional inverse to sigex.zeta2par
	#		bounds: gives bounds for rho and omega, cycle parameters in zeta
	#			rho lies in (bounds[1],bounds[2])
	#			omega lies in (bounds[3],bounds[4])
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Inputs:
	#		mdlPar: see background.  This is the portion of param
	#			corresponding to mdlType, cited as param[[3]]
	#		mdlType: this is a component of mdl (the specified sigex model),
	#			 cited as mdl[[2]]
	#	Outputs:
	#		zeta: see background.
	#	Requires: var.par2pre
	#
	####################################################################

phi2psi <- function(phi)
{
	p <- length(phi)
	pacfs <- phi[p]
	if(p > 1)
	{
		phi <- as.vector(phi[-p])
		for(j in p:2)
		{
			A.mat <- diag(j-1) - pacfs[1]*diag(j-1)[,(j-1):1,drop=FALSE]
			phi <- solve(A.mat,phi)
			pacfs <- c(phi[j-1],pacfs)
			phi <- phi[-(j-1)]
		}
	}
	psi <- log(1+pacfs) - log(1-pacfs)
	return(psi)
}

	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]
	mdlBounds <- mdlType[[3]]

	#############################
	## get zeta for the component

	# ARMA
	if(mdlClass %in% c("arma","arma.stab"))
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		zeta.ar <- NULL
		zeta.ma <- NULL
		if(p.order > 0)
		{
			ar.coef <- mdlPar[1:p.order]
			zeta.ar <- phi2psi(ar.coef)
		}
		if(q.order > 0)
		{
			ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
			zeta.ma <- phi2psi(-1*ma.coef)
		}
		zeta <- c(zeta.ar,zeta.ma)
	}

	# SARMA
	if(mdlClass %in% c("sarma","sarma.stab"))
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		zeta.ar <- NULL
		zeta.ma <- NULL
		zeta.ars <- NULL
		zeta.mas <- NULL
		if(p.order > 0)
		{
			ar.coef <- mdlPar[1:p.order]
			zeta.ar <- phi2psi(ar.coef)
		}
		if(q.order > 0)
		{
			ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
			zeta.ma <- phi2psi(ma.coef)
		}
		if(ps.order > 0)
		{
			ars.coef <- mdlPar[(p.order+q.order+1):(p.order+q.order+ps.order)]
			zeta.ars <- phi2psi(ars.coef)
		}
		if(qs.order > 0)
		{
			mas.coef <- mdlPar[(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order)]
			zeta.mas <- phi2psi(mas.coef)
		}
		zeta <- c(zeta.ar,zeta.ma,zeta.ars,zeta.mas)
	}

	# VARMA
	if(mdlClass %in% c("varma"))
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		zeta.ar <- NULL
		zeta.ma <- NULL
		if(p.order > 0)
		{
			ar.coef <- mdlPar[,,1:p.order,drop=FALSE]
			zeta.ar <- var.par2pre(ar.coef)
		}
		if(q.order > 0)
		{
			ma.coef <- mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE]
			zeta.ma <- var.par2pre(-1*ma.coef)
		}
		zeta <- c(zeta.ar,zeta.ma)
	}

	# SVARMA
	if(mdlClass %in% c("svarma"))
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		zeta.ar <- NULL
		zeta.ma <- NULL
		zeta.ars <- NULL
		zeta.mas <- NULL
		if(p.order > 0)
		{
			ar.coef <- mdlPar[,,1:p.order,drop=FALSE]
			zeta.ar <- var.par2pre(ar.coef)
		}
		if(q.order > 0)
		{
			ma.coef <- mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE]
			zeta.ma <- var.par2pre(ma.coef)
		}
		if(ps.order > 0)
		{
			ars.coef <- mdlPar[,,(p.order+q.order+1):(p.order+q.order+ps.order),drop=FALSE]
			zeta.ars <- var.par2pre(ars.coef)
		}
		if(qs.order > 0)
		{
			mas.coef <- mdlPar[,,(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order),drop=FALSE]
			zeta.mas <- var.par2pre(mas.coef)
		}
		zeta <- c(zeta.ar,zeta.ma,zeta.ars,zeta.mas)
	}

	# cycles
	if(mdlClass %in% c("bw","bw.stab","bal","bal.stab"))
	{
		low.rho <- mdlBounds[1]
		upp.rho <- mdlBounds[2]
		low.omega <- mdlBounds[3]
		upp.omega <- mdlBounds[4]

		x.rho <- (mdlPar[1] - low.rho)/(upp.rho - low.rho)
		rho <- log(x.rho) - log(1 - x.rho)
		x.omega <- (mdlPar[2] - low.omega)/(upp.omega - low.omega)
		omega <- log(x.omega) - log(1 - x.omega)
		zeta <- c(rho,omega)
	}

	# Damped trend
	if(mdlClass %in% c("damped"))
	{
		low <- mdlBounds[1]
		upp <- mdlBounds[2]
		phi <- (mdlPar[1] - low)/(upp - low)
		ar.coef <- log(phi) - log(1 - phi)
		zeta <- ar.coef
	}

	return(zeta)
}
