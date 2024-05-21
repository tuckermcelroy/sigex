#' Compute the autocovariance function of a differenced latent component
#'
#' Background:
#'		A sigex model consists of process x = sum y, for
#'	stochastic components y.  Each component process y_t
#'		is either stationary or is reduced to stationarity by
#'		application of a differencing polynomial delta(B), i.e.
#'			w_t = delta(B) y_t   is stationary.
#'		We have a model for each w_t process, and can compute its
#'		autocovariance function (acf), and denote its autocovariance
#'		generating function (acgf) via gamma_w (B).
#'			Sometimes we may over-difference,
#'		which means applying a differencing polynomial eta(B) that
#'		contains delta(B) as a factor: eta(B) = delta(B)*nu(B).
#'		Then  eta(B) y_t = nu(B) w_t, and the corresponding
#'		acgf is   nu(B) * nu(B^{-1}) * gamma_w (B).
#'
#'		Notes: this function computes the over-differenced acgf,
#'		it is presumed that the given eta(B) contains the needed delta(B)
#'		for that particular component.
#' Conventions: ARMA and VAR models use minus convention for (V)AR polynomials,
#'   and additive convention for (V)MA polynomials.  SARMA and SVARMA use
#'   minus convention for all polynomials.
#'
#' @param L.par Unit lower triangular matrix in GCD of the component's
#'			white noise covariance matrix.
#' @param D.par Vector of logged entries of diagonal matrix in GCD
#'			of the component's white noise covariance matrix.
#' @param mdl The specified sigex model, a list object
#' @param comp Index of the latent component
#' @param mdlPar This is the portion of param
#'			corresponding to mdl[[2]], cited as param[[3]]
#' @param	delta Differencing polynomial (corresponds to eta(B) in Background)
#'			written in format c(delta0,delta1,...,deltad)
#' @param	maxlag Number of autocovariances required
#' @param freqdom A flag, indicating whether frequency domain acf routine should be used.
#'
#' @return 	x.acf: matrix of dimension N x N*maxlag, consisting of autocovariance
#'			matrices stacked horizontally
#' @export
#'

sigex.acf <-
  function(L.par,
           D.par,
           mdl,
           comp,
           mdlPar,
           delta,
           maxlag,
           freqdom = FALSE
  )
{

	##########################################################################
	#
	#	sigex.acf
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
	#	Purpose: compute the autocovariance function of a differenced latent component
	#	Background:
	#		A sigex model consists of process x = sum y, for
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, and can compute its
	#		autocovariance function (acf), and denote its autocovariance
	#		generating function (acgf) via gamma_w (B).
	#			Sometimes we may over-difference,
	#		which means applying a differencing polynomial eta(B) that
	#		contains delta(B) as a factor: eta(B) = delta(B)*nu(B).
	#		Then  eta(B) y_t = nu(B) w_t, and the corresponding
	#		acgf is   nu(B) * nu(B^{-1}) * gamma_w (B).
	#	Notes: this function computes the over-differenced acgf,
	#		it is presumed that the given eta(B) contains the needed delta(B)
	#		for that particular component.
  # Conventions: ARMA and VAR models use minus convention for (V)AR polynomials,
  #   and additive convention for (V)MA polynomials.  SARMA and SVARMA use
  #   minus convention for all polynomials.
	#	Inputs:
	#		L.par: unit lower triangular matrix in GCD of the component's
	#			white noise covariance matrix.  (Cf. sigex.param2gcd background)
	#		D.par: vector of logged entries of diagonal matrix in GCD
	#			of the component's white noise covariance matrix.
	#			(Cf. sigex.param2gcd background)
	#		mdl: the specified sigex model, a list object
	#		comp: index of the latent component
	#		mdlPar: see background to sigex.par2zeta.  This is the portion of param
	#			corresponding to mdl[[2]], cited as param[[3]]
	#		delta: differencing polynomial (corresponds to eta(B) in Background)
	#			written in format c(delta0,delta1,...,deltad)
 	#		maxlag: number of autocovariances required
  #   freqdom: a flag, indicating whether frequency domain acf routine should be used;
  #     for now the default is FALSE, but this is set to TRUE for SARMA and SVARMA models.
	#	Outputs:
	#		x.acf: matrix of dimension N x N*maxlag, consisting of autocovariance
	#			matrices stacked horizontally, i.e.
	#			x.acf = [ gamma(0), gamma(1), ..., gamma(maxlag-1)]
	#	Requires: polymult, polysum, polymulMat, ARMAauto, VARMAauto, specFact,
	#		specFactmvar, sigex.getcycle, sigex.canonize, ubgenerator
	#
	####################################################################

	mdlType <- mdl[[2]][[comp]]
	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]
	mdlBounds <- mdlType[[3]]
	d.delta <- length(delta)
	xi.mat <- L.par %*% diag(exp(D.par),nrow=length(D.par)) %*% t(L.par)
  N <- dim(L.par)[1]

	##################################
	## get acf of stationary component

	# ARMA model
	if(mdlClass == "arma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0)
		{
		  for(j in 1:p.order)
		  {
		    ar.coef <- cbind(ar.coef,diag(mdlPar[,j],nrow=N))
		  }
		}
		if(q.order > 0)
		{
		  for(j in 1:q.order)
		  {
		    ma.coef <- cbind(ma.coef,diag(mdlPar[,j+p.order],nrow=N))
		  }
		}
		ma.array <- array(cbind(diag(N),ma.coef),c(N,N,q.order+1))
		delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
		madiff.array <- polymulMat(delta.array,ma.array)
		ma.coef <- matrix(madiff.array[,,-1],nrow=N)
		if(freqdom)
		{
		  psi.acf <- auto_VARMA(cbind(ar.coef,ma.coef,xi.mat),
		                        p.order,q.order+d.delta-1,0,0,
		                        1,2000,maxlag)[,,1:maxlag,drop=FALSE]
		} else
		{
		  psi.acf <- VARMA_auto(cbind(ar.coef,
		                              matrix(madiff.array[,,-1],nrow=N),
		                              xi.mat),p.order,q.order+d.delta-1,
		                              maxlag)[,,1:maxlag,drop=FALSE]
		}
		x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)
	}

	# Stabilized ARMA model
	if(mdlClass == "arma.stab")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(c(1,-1*ar.coef),canon.delta)
		ma.stab <- sigex.canonize(ma.coef,-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = ar.coef, ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# SARMA model
	if(mdlClass == "sarma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		s.div <- floor(s.period)
		s.frac <- s.period - s.div

		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		ars.coef.stretch <- NULL
		mas.coef.stretch <- NULL
		if(p.order > 0)
		{
		  for(j in 1:p.order)
		  {
		    ar.coef <- cbind(ar.coef,diag(mdlPar[,j],nrow=N))
		  }
		}
		if(q.order > 0)
		{
		  for(j in 1:q.order)
		  {
		    ma.coef <- cbind(ma.coef,diag(mdlPar[,j+p.order],nrow=N))
		  }
		}
		if(s.frac==0)
		{

  		stretch <- c(rep(0,s.period-1),1)
		  if(ps.order > 0)
		  {
  		  for(j in 1:ps.order)
		    {
  			  ars.coef.stretch <- cbind(ars.coef.stretch,t(stretch) %x% diag(mdlPar[,j+p.order+q.order],nrow=N))
			    ars.coef <- cbind(ars.coef,diag(mdlPar[,j+p.order+q.order],nrow=N))
		    }
		  }
	    if(qs.order > 0)
  	  {
	      for(j in 1:qs.order)
	      {
  		    mas.coef.stretch <- cbind(mas.coef.stretch,t(stretch) %x% diag(mdlPar[,j+p.order+q.order+ps.order],nrow=N))
		      mas.coef <- cbind(mas.coef,diag(mdlPar[,j+p.order+q.order+ps.order],nrow=N))
	      }
	    }

		} else  # s.frac > 0
		{
		  freqdom <- FALSE  # we can only compute acf by time domain method
		  if(ps.order > 0) # then ps.order = 1 for this model
		  {
		    for(k in 1:N)
		    {
		      rho.s <- mdlPar[k,1+p.order+q.order]
		      rho.s <- rho.s^(1/s.period)
		      trunc.len <- floor((s.period-1)/2)
		      sar.op <- ubgenerator(s.period,trunc.len,1000,rho.s)
		      sar.op <- polymult(sar.op,c(1,-1*rho.s))
		      if(s.div %% 2 == 0){
		        sar.op = polymult(sar.op, c(1, 1 * rho.s))
		      }
		      ars.coef.stretch <- rbind(ars.coef.stretch,-1*sar.op[-1])
		    }
		  }
		  if(qs.order > 0) # then qs.order = 1 for this model
		  {
		    for(k in 1:N)
		    {
		      rho.s <- mdlPar[k,1+p.order+q.order+ps.order]
		      rho.s <- rho.s^(1/s.period)
		      trunc.len <- floor((s.period-1)/2)
		      sma.op <- ubgenerator(s.period,trunc.len,1000,rho.s)
		      sma.op <- polymult(sma.op,c(1,-1*rho.s))
		      if(s.div %% 2 == 0){
		        sma.op = polymult(sma.op, c(1, 1 * rho.s))
		      }
		      mas.coef.stretch <- rbind(mas.coef.stretch,-1*sma.op[-1])
		    }
		  }


		}

	  delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
    madiff.array <- polymulMat(delta.array,array(cbind(diag(N),-1*ma.coef),c(N,N,q.order+1)))
    if(freqdom)
    {
      ma.coef <- matrix(madiff.array[,,-1],nrow=N)
      psi.acf <- auto_VARMA(cbind(ar.coef,ma.coef,ars.coef,-1*mas.coef,xi.mat),
		                            p.order,q.order+d.delta-1,ps.order,qs.order,
		                            s.period,5000,maxlag)[,,1:maxlag,drop=FALSE]
    } else
    {
      ar.array <- array(cbind(diag(N),-1*ar.coef),c(N,N,p.order+1))
      ars.array <- cbind(diag(N),-1*ars.coef.stretch)
      ars.array <- array(ars.array,c(N,N,dim(ars.array)[2]))
      ar.poly <- polymulMat(ar.array,ars.array)
      ar.coef <- matrix(-1*ar.poly[,,-1],nrow=N)
      mas.array <- cbind(diag(N),-1*mas.coef.stretch)
      mas.array <- array(mas.array,c(N,N,dim(mas.array)[2]))
      ma.poly <- polymulMat(madiff.array,mas.array)
      ma.coef <- matrix(ma.poly[,,-1],nrow=N)
      psi.acf <- VARMA_auto(cbind(ar.coef,ma.coef,xi.mat),
                            dim(ar.poly)[3]-1,
                            dim(ma.poly)[3]-1,
                            maxlag)[,,1:maxlag,drop=FALSE]
    }
		x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)
	}

  # Fractional SARMA of JDemetra+; need ps, qs <= 1
  if(mdlClass == "sarmaf")
  {
    p.order <- mdlOrder[1]
    q.order <- mdlOrder[2]
    ps.order <- mdlOrder[3]
    qs.order <- mdlOrder[4]
    s.period <- mdlOrder[5]
    s.div <- floor(s.period)
    s.frac <- s.period - s.div

    ar.coef <- NULL
    ma.coef <- NULL
    ars.coef <- NULL
    mas.coef <- NULL
    ars.coef.stretch <- NULL
    mas.coef.stretch <- NULL
    freqdom <- FALSE  # we can only compute acf by time domain method
    if(p.order > 0)
    {
      for(j in 1:p.order)
      {
        ar.coef <- cbind(ar.coef,diag(mdlPar[,j],nrow=N))
      }
    }
    if(q.order > 0)
    {
      for(j in 1:q.order)
      {
        ma.coef <- cbind(ma.coef,diag(mdlPar[,j+p.order],nrow=N))
      }
    }
    if(ps.order > 0) # then ps.order = 1 for this model
    {
      for(k in 1:N)
      {
        rho.s <- mdlPar[k,1+p.order+q.order]
        sar.op <- c(1,rep(0,s.div-1),(s.frac-1)*rho.s,-1*s.frac*rho.s)
        ars.coef.stretch <- rbind(ars.coef.stretch,-1*sar.op[-1])
      }
    }
    if(qs.order > 0) # then qs.order = 1 for this model
    {
      for(k in 1:N)
      {
        rho.s <- mdlPar[k,1+p.order+q.order+ps.order]
        sma.op <- c(1,rep(0,s.div-1),(s.frac-1)*rho.s,-1*s.frac*rho.s)
        mas.coef.stretch <- rbind(mas.coef.stretch,-1*sma.op[-1])
      }
    }

    delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
    madiff.array <- polymulMat(delta.array,array(cbind(diag(N),-1*ma.coef),c(N,N,q.order+1)))

      ar.array <- array(cbind(diag(N),-1*ar.coef),c(N,N,p.order+1))
      ars.array <- cbind(diag(N),-1*ars.coef.stretch)
      ars.array <- array(ars.array,c(N,N,dim(ars.array)[2]))
      ar.poly <- polymulMat(ar.array,ars.array)
      ar.coef <- matrix(-1*ar.poly[,,-1],nrow=N)
      mas.array <- cbind(diag(N),-1*mas.coef.stretch)
      mas.array <- array(mas.array,c(N,N,dim(mas.array)[2]))
      ma.poly <- polymulMat(madiff.array,mas.array)
      ma.coef <- matrix(ma.poly[,,-1],nrow=N)
      psi.acf <- VARMA_auto(cbind(ar.coef,ma.coef,xi.mat),
                            dim(ar.poly)[3]-1,
                            dim(ma.poly)[3]-1,
                            maxlag)[,,1:maxlag,drop=FALSE]

    x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)

  }


	# Stabilized SARMA model
	if(mdlClass == "sarma.stab")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		stretch <- c(rep(0,s.period-1),1)
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		if(ps.order > 0)
		{
			ars.coef <- mdlPar[(p.order+q.order+1):(p.order+q.order+ps.order)]
			ars.coef <- ars.coef %x% stretch
		}
		if(qs.order > 0)
		{
			mas.coef <- mdlPar[(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order)]
			mas.coef <- mas.coef %x% stretch
		}
		ar.poly <- polymult(c(1,-1*ar.coef),c(1,-1*ars.coef))
		ma.poly <- polymult(c(1,-1*ma.coef),c(1,-1*mas.coef))
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# VARMA model
	if(mdlClass == "varma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- matrix(mdlPar[,,1:p.order,drop=FALSE],nrow=N)
		if(q.order > 0) ma.coef <- matrix(mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE],nrow=N)
		delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
		madiff.array <- polymulMat(delta.array,array(cbind(diag(N),ma.coef),c(N,N,q.order+1)))
		ma.coef <- matrix(madiff.array[,,-1],nrow=N)
		if(freqdom)
		{
		  psi.acf <- auto_VARMA(cbind(ar.coef,ma.coef,xi.mat),
		                        p.order,q.order+d.delta-1,0,0,
		                        1,2000,maxlag)[,,1:maxlag,drop=FALSE]
		} else
		{
		  psi.acf <- VARMA_auto(cbind(ar.coef,ma.coef,xi.mat),
		                        p.order,q.order+d.delta-1,
		                        maxlag)[,,1:maxlag,drop=FALSE]
		}
		x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)
	}

	# SVARMA model
	if(mdlClass == "svarma")
	{
#	  freqdom <- TRUE
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		stretch <- c(rep(0,s.period-1),1)
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		ars.coef.stretch <- NULL
		mas.coef.stretch <- NULL
		if(p.order > 0)
		{
			ar.coef <- matrix(mdlPar[,,1:p.order,drop=FALSE],nrow=N)
		}
		if(q.order > 0)
		{
			ma.coef <- matrix(mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE],nrow=N)
		}
		if(ps.order > 0)
		{
			ars.coef <- matrix(mdlPar[,,(p.order+q.order+1):(p.order+q.order+ps.order),drop=FALSE],nrow=N)
			for(j in 1:ps.order)
			{
			  ars.coef.stretch <- cbind(ars.coef.stretch,t(stretch) %x% mdlPar[,,j+p.order+q.order])
			}
		}
		if(qs.order > 0)
		{
			mas.coef <- matrix(mdlPar[,,(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order),drop=FALSE],nrow=N)
			for(j in 1:qs.order)
			{
			  mas.coef.stretch <- cbind(mas.coef.stretch,t(stretch) %x% mdlPar[,,j+p.order+q.order+ps.order])
			}
		}
		delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
		madiff.array <- polymulMat(delta.array,array(cbind(diag(N),-1*ma.coef),c(N,N,q.order+1)))
		if(freqdom)
		{
		  ma.coef <- matrix(madiff.array[,,-1],nrow=N)
		  psi.acf <- auto_VARMA(cbind(ar.coef,ma.coef,ars.coef,-1*mas.coef,xi.mat),
		                      p.order,q.order+d.delta-1,ps.order,qs.order,
		                      s.period,5000,maxlag)[,,1:maxlag,drop=FALSE]
		} else
		{
		  ar.array <- array(cbind(diag(N),-1*ar.coef),c(N,N,p.order+1))
		  ars.array <- array(cbind(diag(N),-1*ars.coef.stretch),c(N,N,s.period*ps.order+1))
		  ar.poly <- polymulMat(ar.array,ars.array)
		  ar.coef <- matrix(-1*ar.poly[,,-1],nrow=N)
		  mas.array <- array(cbind(diag(N),-1*mas.coef.stretch),c(N,N,s.period*qs.order+1))
		  ma.poly <- polymulMat(madiff.array,mas.array)
		  ma.coef <- matrix(ma.poly[,,-1],nrow=N)
		  psi.acf <- VARMA_auto(cbind(ar.coef,ma.coef,xi.mat),
		                        p.order+s.period*ps.order,
		                        q.order+d.delta-1+s.period*qs.order,
		                        maxlag)[,,1:maxlag,drop=FALSE]
		}
		x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)
	}

	# Butterworth cycle
	if(mdlClass == "bw")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		ma.poly <- out[[2]]
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1],ma = polymult(ma.poly,delta)[-1],
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized Butterworth cycle
	if(mdlClass == "bw.stab")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		ma.poly <- out[[2]]
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Balanced cycle
	if(mdlClass == "bal")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}
		ma.acf <- c(rev(ma.acf),ma.acf[-1]) + 1e-10
		ma.poly <- Re(specFact(ma.acf))
		ma.scale <- ma.poly[1]^2
		ma.poly <- ma.poly/ma.poly[1]
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1],ma = polymult(ma.poly,delta)[-1],
			lag.max=maxlag)[1:maxlag]
	 	psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized Balanced cycle
	if(mdlClass == "bal.stab")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}
		ma.acf <- c(rev(ma.acf),ma.acf[-1]) + 1e-10
		ma.poly <- Re(specFact(ma.acf))
		ma.scale <- ma.poly[1]^2
		ma.poly <- ma.poly/ma.poly[1]
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.scale*ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Damped Trend model
	if(mdlClass == "damped")
	{
		p.order <- mdlOrder[1]
		ar.coef <- mdlPar[1]
		ar.poly <- 1
		for(k in 1:p.order)
		{
			ar.poly <- polymult(ar.poly,c(1,-1*ar.coef))
		}
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = delta[-1],lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}

  # ARMA Copula
  if(mdlClass == "armacopula")
  {
    p.order <- mdlOrder[1,]
    q.order <- mdlOrder[2,]
    p.max <- max(p.order)
    q.max <- max(q.order)
    ar.coef <- NULL
    ma.coef <- NULL
    if(p.max > 0)
    {
      for(j in 1:p.max)
      {
        ar.coef <- cbind(ar.coef,diag(mdlPar[,j],nrow=N))
      }
    }
    if(q.max > 0)
    {
      for(j in 1:q.max)
      {
        ma.coef <- cbind(ma.coef,diag(mdlPar[,j+p.max],nrow=N))
      }
    }
    ma.array <- array(cbind(diag(N),ma.coef),c(N,N,q.max+1))
    delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
    madiff.array <- polymulMat(delta.array,ma.array)
    ma.coef <- matrix(madiff.array[,,-1],nrow=N)
    if(freqdom)
    {
      psi.acf <- auto_VARMA(cbind(ar.coef,ma.coef,xi.mat),
                            p.order,q.order+d.delta-1,0,0,
                            1,2000,maxlag)[,,1:maxlag,drop=FALSE]
    } else
    {
      psi.acf <- VARMA_auto(cbind(ar.coef,
                                  matrix(madiff.array[,,-1],nrow=N),
                                  xi.mat),p.order,q.order+d.delta-1,
                            maxlag)[,,1:maxlag,drop=FALSE]
    }
    x.acf <- matrix(aperm(psi.acf,c(1,3,2)),ncol=N)
  }

	return(x.acf)
}

