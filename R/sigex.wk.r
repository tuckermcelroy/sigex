sigex.wk <- function(data.ts,param,mdl,sigcomps,target,plotit=TRUE,grid,len)
{

	##########################################################################
	#
	#	sigex.wk
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
	#	Purpose: computes signal extraction filter coefficients and MSE
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
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Notes: take grid >> len, else numerical issues arise
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#		sigcomps: indices of the latent components composing the signal S_t
	#   target: array of dimension c(N,N,p+1) for degree p matrix polynomial
  #     phi (B) such that target signal is phi (B) S_t.
  #     normal usage is target <- array(diag(N),c(N,N,1))
  #		plotit: binary flag to indicate whether plots should be generated;
	#			only gives output if N <= 3
	#		grid: desired number of frequencies for spectrum calculations
	#		len: max index of the filter coefficients
	#	Outputs:
	#		list object of psi.wk and mse.wksig
	#		psi.wk: array of dimension c(N,N,2*len+1), 
	#			of real number entries, where the third dimension
	#			indexes coefficients from -len,...,0,...,len 
	#		mse.wksig: N x N matrix of signal extraction mean squared error
	#	Requires: sigex.frf, sigex.wkmse
	#
	####################################################################
 
quad <- function(z)
{
	# quadrature of z 
	len <- length(z)
	out <- (z[1]+z[len])/2 + sum(z[2:(len-1)])
	return(out/len)
}

	x <- t(data.ts)
	N <- dim(x)[1]
	lambda <- pi*seq(0,grid)/grid
	
	f.phi <- t(rep(1,(grid+1))) %x% target[,,1]  
	if(dim(target)[3] > 1) {
	  for(i in 2:dim(target)[3])
	  {
	    f.phi <- f.phi + t(exp(-1i*lambda*(i-1))) %x% target[,,i]  
	  } }
	f.phi <- array(f.phi,c(N,N,(grid+1)))
	
	### wk filter 
	frf.wk <- sigex.frf(data.ts,param,mdl,sigcomps,grid)
	frf.new <- do.call(cbind,lapply(seq(1,grid+1),function(i) f.phi[,,i] %*% frf.wk[,,i] ))
	frf.wk <- array(frf.new,c(N,N,grid+1))
	psi.wk <- Re(apply(frf.wk,c(1,2),quad))
	for(i in 1:len)
	{
		integrand <- t(exp(1i*lambda*i) %x% matrix(1,N,N)) * 
			matrix(frf.wk,nrow=N)
		integrand <- array(integrand,c(N,N,(grid+1)))
		next.coeff <- Re(apply(integrand,c(1,2),quad))
		psi.wk <- cbind(psi.wk,next.coeff)
		integrand <- t(exp(-1i*lambda*i) %x% matrix(1,N,N)) * 
			matrix(frf.wk,nrow=N)
		integrand <- array(integrand,c(N,N,(grid+1)))
		prev.coeff <- Re(apply(integrand,c(1,2),quad))
		psi.wk <- cbind(prev.coeff,psi.wk)
	}
	psi.wk <- array(psi.wk,c(N,N,(2*len+1)))

	### wk MSE
	frf.wksig <- sigex.wkmse(data.ts,param,mdl,sigcomps,grid)
	frf.new <- do.call(cbind,lapply(seq(1,grid+1),function(i) f.phi[,,i] %*% frf.wksig[,,i] %*% Conj(t(f.phi[,,i]))))
	frf.wksig <- array(frf.new,c(N,N,grid+1))
	mse.wksig <- Re(apply(frf.wksig,c(1,2),quad))

	if(plotit && (N <= 4)) {
	par(mfrow = c(N,N))
	for(i in 1:N) 
	{
		for(j in 1:N)
		{
			plot(ts(psi.wk[i,j,],start=-len,frequency=1),
				xlab="Index",ylab="")
		}	
	} }

	return(list(psi.wk,mse.wksig))
}


