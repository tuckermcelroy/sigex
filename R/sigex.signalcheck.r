sigex.signalcheck <- function(signal,param,mdl,sigcomps,lagall)
{

	##########################################################################
	#
	#	sigex.signalcheck
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
	#	Purpose: computes model-based signal extraction diagnostics
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
	#		Signal extraction diagnostics are computed, which are based
	#		on autocorrelations of differenced estimated signals, compared
	#		to their expected variability.  Blakely and McElroy (2017 Econometric Reviews)
	#		discusses the theory in generality; this version does
	#		not take parameter uncertainty into account.
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		signal: a T x N matrix ts object of the extracted signal
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#		sigcomps: indices of the latent components composing the signal
	#		lagall: each signal diagnostic has a lag argument.  They are all
	#			computed, from lags one through lagall
	#	Outputs:
	#		diagnostics: complex matrix with lagall rows and N columns, one
	#			column for each series.  The real part is the diagnostic
	#			value, given by sample autocorrelation at various lags of
	#			the estimated signal.  The imaginary part is the t statistic
	#			of the diagnostic: statistic minus null mean, divided by
	#			standard error.
	#	Requires: sigex.delta, sigex.acf
	#
	####################################################################
 
	x <- t(signal)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# compute differencing polynomials and matrices
	allcomps <- seq(1,length(mdl[[3]]))
	noisecomps <- allcomps[!allcomps %in% sigcomps]
	delta.signal <- sigex.delta(mdl,noisecomps)
	TdiffSig <- T - length(delta.signal) + 1

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	acfsignal.mat <- matrix(0,nrow=N*TdiffSig,ncol=N)
	for(i in sigcomps)
	{
		L.par[[i]] <- param[[1]][[i]]
		D.par[[i]] <- param[[2]][[i]]
		mdlType <- mdl[[2]][i]
		delta <- sigex.delta(mdl,c(noisecomps,i))
		acfsignal.mat <- acfsignal.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,param[[3]][[i]],delta,TdiffSig)		
	}
	signal.acf <- array(acfsignal.mat,dim=c(N,TdiffSig,N))

	dsig <- filter(signal,delta.signal,method="convolution",sides=1)[length(delta.signal):T,]
 	acf.sample <- acf(dsig,type="correlation",lag.max=TdiffSig,plot=FALSE)$acf
	
	L <- TdiffSig-1
	diagnostics <- NULL
	for(i in 1:N)
	{
		auto.corr <- signal.acf[i,,i]/signal.acf[i,1,i]
		sig.corrs <- auto.corr[2:(lagall+1)]
		est.corrs <- acf.sample[2:(lagall+1),i,i]
		auto.corr <- c(rev(auto.corr),auto.corr[-1])
		all.w <- NULL

		for(h in 1:lagall) 
		{
			new.w <- sum((auto.corr[(L+1+1+h):(L+1+L)] + auto.corr[(L+1+1-h):(L+1+L-2*h)] -
				2*auto.corr[L+1+h]*auto.corr[(L+1+1):(L+1+L-h)])^2)
			all.w <- c(all.w,new.w)
		}
		ses <- sqrt(abs(all.w)/TdiffSig)	
		diagnostics <- cbind(diagnostics,est.corrs + 1i*(est.corrs - sig.corrs)/ses)
	}

	return(diagnostics)
}
