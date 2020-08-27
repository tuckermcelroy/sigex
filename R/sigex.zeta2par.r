sigex.zeta2par <- function(zeta,mdlType)
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
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Notes: this is a functional inverse to sigex.par2zeta
	#		bounds: gives bounds for rho and omega, cycle parameters in zeta
	#			rho lies in (bounds[1],bounds[2])
	#			omega lies in (bounds[3],bounds[4])
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Inputs:
	#		zeta: see background.  This is the portion of psi
	#			corresponding t.s. models, such as cycles
	#		mdlType: this is a component of mdl (the specified sigex model),
	#			 cited as mdl[[2]]
	#	Outputs:
	#		zeta.par: see background.  This is a portion of the full
	#			param list, corresponding to param[[3]]
	#	Requires: sigex.param2gcd, var.pre2par
	#
	####################################################################

psi2phi <- function(psi)
{
	p <- length(psi)
	pacfs <- (exp(psi)-1)/(exp(psi)+1)
	if(p==1)
	{
		phi <- pacfs
	} else
	{
		phi <- as.vector(pacfs[1])
		for(j in 2:p)
		{
			A.mat <- diag(j-1) - pacfs[j]*diag(j-1)[,(j-1):1,drop=FALSE]
			phi <- A.mat %*% phi
			phi <- c(phi,pacfs[j])
		}
	}
	return(phi)
}

	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]
	mdlBounds <- mdlType[[3]]

	##############################
	## get param for the component

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
			zeta.ar <- zeta[1:p.order]
			ar.coef <- psi2phi(zeta.ar)
		}
		if(q.order > 0)
		{
			zeta.ma <- zeta[(p.order+1):(p.order+q.order)]
			ma.coef <- psi2phi(-1*zeta.ma)
		}
		zeta.par <- c(ar.coef,ma.coef)
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
			zeta.ar <- zeta[1:p.order]
			ar.coef <- psi2phi(zeta.ar)
		}
		if(q.order > 0)
		{
			zeta.ma <- zeta[(p.order+1):(p.order+q.order)]
			ma.coef <- psi2phi(zeta.ma)
		}
		if(ps.order > 0)
		{
			zeta.ars <- zeta[(p.order+q.order+1):(p.order+q.order+ps.order)]
			ars.coef <- psi2phi(zeta.ars)
		}
		if(qs.order > 0)
		{
			zeta.mas <- zeta[(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order)]
			mas.coef <- psi2phi(zeta.mas)
		}
		zeta.par <- c(ar.coef,ma.coef,ars.coef,mas.coef)
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
		zeta.par <- NULL
		if(p.order > 0)
		{
			zeta.ar <- zeta[1:(p.order*N^2)]
			ar.coef <- var.pre2par(zeta.ar,p.order,N,FALSE)
			zeta.par <- matrix(ar.coef,nrow=N)
		}
		if(q.order > 0)
		{
			zeta.ma <- zeta[(p.order*N^2 +1):(p.order*N^2 + q.order*N^2)]
			ma.coef <- -1*var.pre2par(zeta.ma,q.order,N,FALSE)
			zeta.par <- cbind(zeta.par,matrix(ma.coef,nrow=N))
		}
		zeta.par <- array(zeta.par,c(N,N,p.order+q.order))
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
		zeta.par <- NULL
		if(p.order > 0)
		{
			zeta.ar <- zeta[1:(p.order*N^2)]
			ar.coef <- var.pre2par(zeta.ar,p.order,N,FALSE)
			zeta.par <- matrix(ar.coef,nrow=N)
		}
		if(q.order > 0)
		{
			zeta.ma <- zeta[(p.order*N^2 +1):(p.order*N^2 + q.order*N^2)]
			ma.coef <- var.pre2par(zeta.ma,q.order,N,FALSE)
			zeta.par <- cbind(zeta.par,matrix(ma.coef,nrow=N))
		}
		if(ps.order > 0)
		{
			zeta.ars <- zeta[(p.order*N^2+q.order*N^2+1):(p.order*N^2+q.order*N^2+ps.order*N^2)]
			ars.coef <- var.pre2par(zeta.ars,ps.order,N,FALSE)
			zeta.par <- cbind(zeta.par,matrix(ars.coef,nrow=N))
		}
		if(qs.order > 0)
		{
			zeta.mas <- zeta[(p.order*N^2+q.order*N^2+ps.order*N^2+1):(p.order*N^2+q.order*N^2+ps.order*N^2+qs.order*N^2)]
			mas.coef <- var.pre2par(zeta.mas,qs.order,N,FALSE)
			zeta.par <- cbind(zeta.par,matrix(mas.coef,nrow=N))
		}
		zeta.par <- array(zeta.par,c(N,N,p.order+q.order+ps.order+qs.order))
	}

	# cycles
	if(mdlClass %in% c("bw","bw.stab","bal","bal.stab"))
	{
		low.rho <- mdlBounds[1]
		upp.rho <- mdlBounds[2]
		low.omega <- mdlBounds[3]
		upp.omega <- mdlBounds[4]

		rho <- low.rho + (upp.rho - low.rho)*exp(zeta[1])/(1+exp(zeta[1]))
		omega <- low.omega + (upp.omega - low.omega)*exp(zeta[2])/(1+exp(zeta[2]))
		zeta.par <- c(rho,omega)
	}

	# Damped trend
	if(mdlClass %in% c("damped"))
	{
		low <- mdlBounds[1]
		upp <- mdlBounds[2]
		ar.coef <- low + (upp - low)*exp(zeta[1])/(1+exp(zeta[1]))
		zeta.par <- ar.coef
	}

 	return(zeta.par)
}
