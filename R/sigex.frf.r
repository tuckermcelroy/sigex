sigex.frf <- function(data.ts,param,mdl,sigcomps,grid)
{

	##########################################################################
	#
	#	sigex.frf
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
	#	Purpose: computes signal extraction filter frequency response function
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
	#		sigcomps: indices of the latent components composing the signal
	#		grid: desired number of frequencies for spectrum calculations
	#	Outputs:
	#		frf.wk:  array of dimension c(N,N,grid), with complex number entries 
	#	Requires: sigex.spectra, sigex.delta
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	f.all <- t(rep(0,grid+1) %x% diag(N))
	frf.wk <- array(f.all,c(N,N,(grid+1)))
	f.sig <- f.all
	for(i in 1:length(mdl[[3]]))
	{
		L.par <- param[[1]][[i]]
		D.par <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],delta,grid)
		f.all <- f.all + matrix(f.comp,nrow=N)
		if(i %in% sigcomps) f.sig <- f.sig + matrix(f.comp,nrow=N)
	}
	f.all  <- array(f.all,c(N,N,(grid+1)))
	f.sig  <- array(f.sig,c(N,N,(grid+1)))
	for(j in 1:(grid+1))
	{
		flag.zero <- FALSE
		for(k in 1:length(mdl[[3]]))
		{
			delta <- mdl[[3]][[k]]
			if(Mod(sum(delta*exp(-1i*seq(0,length(delta)-1)*(j-1)*pi/grid))) < 10^(-8))
			{
				flag.zero <- TRUE
				g.comp <- 0*diag(N)
				for(i in setdiff(seq(1,length(mdl[[3]])),k))
				{
					L.par <- param[[1]][[i]]
					D.par <- param[[2]][[i]]
					delta <- sigex.delta(mdl,c(i,k))
					f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],delta,grid)[,,j]
					g.comp <- g.comp + matrix(f.comp,nrow=N)
				}
				L.par <- param[[1]][[k]]
				D.par <- param[[2]][[k]]
				h.comp <- t(L.par) %*% solve(g.comp) %*% L.par %*% diag(exp(D.par),nrow=length(mdl[[1]][[k]]))
				h.comp <- L.par %*% diag(exp(D.par),nrow=length(mdl[[1]][[k]])) %*%
						solve(h.comp) %*% t(L.par) %*% solve(g.comp)			
				h.comp <- diag(N) - h.comp
				for(l in sigcomps)
				{
					if(l == k) { sig.comp <- diag(N) - h.comp } else { 
						L.par <- param[[1]][[l]]
						D.par <- param[[2]][[l]]
						delta <- sigex.delta(mdl,c(l,k))
						f.comp <- sigex.spectra(L.par,D.par,mdl,l,param[[3]][[l]],delta,grid)[,,j]		
						sig.comp <- f.comp %*% solve(g.comp) %*% h.comp
					}
					frf.wk[,,j] <- frf.wk[,,j] + sig.comp
				}
			}  
		}
	 	if(!flag.zero) { frf.wk[,,j] <- f.sig[,,j] %*% solve(f.all[,,j]) } 
	}	

	return(frf.wk)
}


