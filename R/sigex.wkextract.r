sigex.wkextract <- function(psi,mdl,data.ts,sigcomps,target,grid,window,horizon,needMSE)
{

	##########################################################################
	#
	#	sigex.wkextract
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

	################# Documentation ############################################
	#
	#	Purpose: computes signal extractions and MSE via WK method
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
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter manifold) 
	#	Notes:
	#		method does midcasts for specified time indices, 
	#		applies WK filter for signal (given by sigcomps),
	#		based on grid for Riemann integration,
	#		truncated to length given by 2*window + 1;
	#		generates signal at all time points t in seq(1-horizon,T+horizon).
	#	Inputs:
	#		psi: see background.  
	#		mdl: the specified sigex model, a list object
  #		data.ts: a T x N matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #		sigcomps: indices of the latent components composing the signal S_t
  #   target: array of dimension c(N,N,p+1) for degree p matrix polynomial
  #     phi (B) such that target signal is phi (B) S_t.
  #     normal usage is target <- array(diag(N),c(N,N,1))
  #		grid: desired number of frequencies for spectrum calculations
	#		window: max index of the WK filter coefficients
	#		horizon: a positive integer indicating how many forecasts and
	#			aftcasts of the signal should be generated
	#		needMSE: a binary flag, set to 1 if you want MSE based on casting error,
	#			or if there are any missing values; else (with value 0) the routine
	#			runs faster and returns only WK portion of MSE.
	#	Outputs:
	#		list object of extract.sig, upp, and low
	#		extract.sig: (T+H) x N matrix of the signal estimates, where H is
	#			twice the length of horizon
	#		upp: as extract.sig, plus twice the standard error
	#		low: as extract.sig, minus twice the standard error
	#	Requires: sigex.psi2par, sigex.wk, sigex.midcast, sigex.cast
	#
	#####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	param <- sigex.psi2par(psi,mdl,data.ts)

	wk.out <- sigex.wk(data.ts,param,mdl,sigcomps,target,FALSE,grid,window)
	wk.filter <- wk.out[[1]]
	wk.mse <- wk.out[[2]]
	wk.len <- 2*window+1
	
	leads.mid <- NULL
	for(k in 1:N) { leads.mid <- union(leads.mid,seq(1,T)[is.na(data.ts)[,k]]) }
	leads.aft <- seq(1-horizon-window,0)
	leads.fore <- seq(T+1,T+horizon+window)
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
		casts.all <- sigex.midcast(psi,mdl,data.ts,horizon+window)
		casts.var <- array(casts.all[[2]],c(N,len.all,N,len.all))
		data.ext <- data.demean
		if(len.mid>0) data.ext[leads.mid,] <- 
			t(casts.all[[1]][,(len.aft+1):(len.aft+len.mid)])
		data.ext <- rbind(t(casts.all[[1]][,1:len.aft,drop=FALSE]),data.ext,
			t(casts.all[[1]][,(len.aft+len.mid+1):len.all,drop=FALSE]))
		t.start <- 1
		for(t in seq(1-horizon,T+horizon))
		{
#			print(t.start)
			range.t <- intersect(seq(t-window,t+window),leads.all)
			len.t <- length(range.t)
			if(len.t==0) { cast.mse[,,t+horizon] <- 0*diag(N) } else {
      wk.coefs <- wk.filter[,,window+1+t-range.t,drop=FALSE]
			cast.mse[,,t+horizon] <- matrix(wk.coefs,c(N,N*len.t)) %*% 
				matrix(casts.var[,seq(0,len.t-1)+t.start,,seq(0,len.t-1)+t.start,drop=FALSE],c(len.t*N,len.t*N)) %*%
				t(matrix(wk.coefs,c(N,N*len.t))) }
			if(is.element(t-window,leads.all)) { t.start <- t.start + 1 }
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
			output.k <- filter(data.ext[,k],wk.filter[j,k,],
				method="convolution",sides=2)
			output.k <- output.k[(window+1):(window+T+2*horizon)]
			output.j <- output.j + output.k
		}
		mse <- wk.mse[j,j] + cast.mse[j,j,]
		extract.sig <- cbind(extract.sig,output.j)
		upp <- cbind(upp,output.j + 2*sqrt(mse))
		low <- cbind(low,output.j - 2*sqrt(mse))
	}	

	return(list(extract.sig,upp,low))
}


