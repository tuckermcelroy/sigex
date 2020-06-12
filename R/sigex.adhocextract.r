sigex.adhocextract <- function(psi,mdl,data.ts,adhoc,shift,horizon,needMSE)
{

  ##########################################################################
  #
  #	sigex.adhocextract
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

  ################# Documentation ############################################
  #
  #	Purpose: computes signal extractions and MSE from an ad hoc filter
  #	Background:
  #		A sigex model consists of process x = sum y, for
  #		stochastic components y.  Each component process y_t
  #		is either stationary or is reduced to stationarity by
  #		application of a differencing polynomial delta(B), i.e.
  #			w_t = delta(B) y_t   is stationary.
  #		We have a model for each w_t process, and can compute its
  #		autocovariance function (acf), and denote its autocovariance
  #		generating function (acgf) via gamma_w (B).
  #		The signal extraction filter for y_t is determined from
  #		this acgf and delta.
  #		param is the name for the model parameters entered into
  #		a list object with a more intuitive structure, whereas
  #		psi refers to a vector of real numbers containing all
  #		hyper-parameters (i.e., reals mapped bijectively to the parameter
  #		manifold) together with imaginary component flagging
  #		whether the hyper-parameter is fixed for purposes of estimation.
  #	Notes:
  #		method does midcasts for specified time indices,
  #		applies ad hoc filter for signal (given by adhoc),
  #		of length L and offset "shift".
  #		generates signal at all time points t in seq(L-shift-horizon,T-shift+horizon).
  #	Inputs:
  #		psi: see background.
  #		mdl: the specified sigex model, a list object
  #		data.ts: a T x N matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   adhoc:  an array N x N x L, where L is length
  #   shift: gives the integer offset for the adhoc filter:
  #     filter coefficients have indices -shift,...,0,...,L-1-shift
  #     set shift = 0 for a causal filter
  #		horizon: a non-negative integer indicating how many forecasts and
  #			aftcasts of the signal should be generated
  #		needMSE: a binary flag, set to 1 if you want MSE based on casting error,
  #			or if there are any missing values; else (with value 0) the routine
  #			runs faster and returns zero for the MSE.
  #	Outputs:
  #		list object of extract.sig, upp, and low
  #		extract.sig: (T+H) x N matrix of the signal estimates, where H is
  #			twice the length of horizon
  #		upp: as extract.sig, plus twice the standard error
  #		low: as extract.sig, minus twice the standard error
  #	Requires: sigex.psi2par, sigex.midcast, sigex.cast
  #
  #####################################################################

  x <- t(data.ts)
  N <- dim(x)[1]
  T <- dim(x)[2]
  param <- sigex.psi2par(psi,mdl,data.ts)
  L <- dim(adhoc)[3]

  # The filter output is sum_{j=-shift}^{L-1-shift} psi_j x_{t-j}.
  # In order to generate output for times t = (1-horizon):(T+horizon),
  # we must forecast times (T+1):(T+horizon+shift)
  # and aftcast times (2-horizon-L+shift):0
  leads.mid <- NULL
  for(k in 1:N) { leads.mid <- union(leads.mid,seq(1,T)[is.na(data.ts)[,k]]) }
  leads.aft <- seq(2-horizon-L+shift,0)
  leads.fore <- seq(T+1,T+horizon+shift)
  leads.all <- union(leads.aft,union(leads.mid,leads.fore))
  len.aft <- length(leads.aft)
  len.fore <- length(leads.fore)
  len.mid <- length(leads.mid)
  len.all <- length(leads.all)

  # remove regression effects from raw data
  beta.par <- param[[4]]
  ind <- 0
  data.demean <- data.ts
  data.demean[is.na(data.demean)] <- 1i
  for(k in 1:N)
  {
    reg.mat <- mdl[[4]][[k]]
    len <- dim(reg.mat)[2]
    data.demean[,k] <- data.demean[,k] - reg.mat %*% beta.par[(ind+1):(ind+len)]
    ind <- ind+len
  }

  cast.mse <- array(0,c(N,N,T+2*horizon))
  if(needMSE)
  {
    # need to get len.aft aftcasts and len.fore forecasts; use sigex.midcast to
    #   get len.aft+len.fore aftcasts and forecasts (each side), and then drop
    #   len.fore aftcasts and len.aft forecasts.
    casts.all <- sigex.midcast(psi,mdl,data.ts,len.aft+len.fore)
    casts.x <- casts.all[[1]][,(len.fore+1):(2*len.fore+len.mid+len.aft),drop=FALSE]
    casts.var <- array(casts.all[[2]],c(N,2*len.aft+len.mid+2*len.fore,N,2*len.aft+len.mid+2*len.fore))
    casts.var <- casts.var[,(len.fore+1):(2*len.fore+len.mid+len.aft),,
                           (len.fore+1):(2*len.fore+len.mid+len.aft),drop=FALSE]
    data.ext <- data.demean
    if(len.mid>0) data.ext[leads.mid,] <- t(casts.x[,(len.aft+1):(len.aft+len.mid)])
    data.ext <- rbind(t(casts.x[,1:len.aft,drop=FALSE]),data.ext,
                      t(casts.x[,(len.aft+len.mid+1):len.all,drop=FALSE]))
    t.start <- 1
    for(t in seq(1-horizon,T+horizon))
    {
      #			print(t.start)
      range.t <- intersect(seq(t-(L-1-shift),t+(shift)),leads.all)
      len.t <- length(range.t)
      if(len.t==0) { cast.mse[,,t+horizon] <- 0*diag(N) } else {
        adhoc.coefs <- adhoc[,,shift+1+t-range.t,drop=FALSE]
        cast.mse[,,t+horizon] <- matrix(adhoc.coefs,c(N,N*len.t)) %*%
          matrix(casts.var[,seq(0,len.t-1)+t.start,,seq(0,len.t-1)+t.start,drop=FALSE],c(len.t*N,len.t*N)) %*%
          t(matrix(adhoc.coefs,c(N,N*len.t))) }
      if(is.element(t-(L-1-shift),leads.all)) { t.start <- t.start + 1 }
    }
  } else
  {
    data.ext <- t(sigex.cast(psi,mdl,data.ts,leads.all))
  }

  extract.sig <- NULL
  upp <- NULL
  low <- NULL
  for(j in 1:N)
  {
    output.j <- rep(0,T+2*horizon)
    for(k in 1:N)
    {
      output.k <- filter(data.ext[,k],adhoc[j,k,],
                         method="convolution",sides=1)
      output.k <- output.k[-seq(1,L-1)]
      output.j <- output.j + output.k
    }
    mse <- pmax(0,cast.mse[j,j,])
    extract.sig <- cbind(extract.sig,output.j)
    upp <- cbind(upp,output.j + 2*sqrt(mse))
    low <- cbind(low,output.j - 2*sqrt(mse))
  }

  return(list(extract.sig,upp,low))
}


