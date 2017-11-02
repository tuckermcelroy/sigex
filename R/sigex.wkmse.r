sigex.wkmse <- function(data,param,mdl,sigcomps,grid)
{

	###############################
	#   sigex.wkmse
	#	by Tucker McElroy	
	#
	#	Computes signal extraction mean square error (bi-infinite sample)
	#		param must be in format yielded by sigex.default
	#		sigcomps provides indices of the desired components
	#		grid is desired resolution of frequencies 
	#			(take to be a large prime, e.g., 997)
	#     Output is in array form of dimension c(N,N,grid), 
	#		of complex number entries, giving signal error frf 
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	f.all <- t(rep(0,grid+1) %x% diag(N))
	frf.wk <- array(f.all,c(N,N,(grid+1)))
	f.sig <- f.all
	f.noise <- f.all
	allcomps <- seq(1,length(mdl[[3]]))
	noisecomps <- allcomps[!allcomps %in% sigcomps]
	for(i in 1:length(mdl[[3]]))
	{
		L.par <- param[[1]][[i]]
		D.par <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],N,delta,grid)
		f.all <- f.all + matrix(f.comp,nrow=N)
		if(i %in% sigcomps)
		{
			delta <- sigex.delta(mdl,c(i,noisecomps))
			f.sigcomp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],N,delta,grid)
			f.sig <- f.sig + matrix(f.sigcomp,nrow=N)
		}
		if(i %in% noisecomps)
		{
			delta <- sigex.delta(mdl,c(i,sigcomps))
			f.noisecomp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],N,delta,grid)
			f.noise <- f.noise + matrix(f.noisecomp,nrow=N)
		}
	}
	f.all  <- array(f.all,c(N,N,(grid+1)))
	f.sig  <- array(f.sig,c(N,N,(grid+1)))
	f.noise <- array(f.noise,c(N,N,(grid+1)))
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
				frf.sig <- 0*diag(N)
				for(i in setdiff(seq(1,length(mdl[[3]])),k))
				{
					L.par <- param[[1]][[i]]
					D.par <- param[[2]][[i]]
					delta <- sigex.delta(mdl,c(i,k))
					f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],N,delta,grid)[,,j]
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
					if(l == k) { sig.comp <- (diag(N) - h.comp) } else { 
						L.par <- param[[1]][[l]]
						D.par <- param[[2]][[l]]
						delta <- sigex.delta(mdl,c(l,k))
						f.comp <- sigex.spectra(L.par,D.par,mdl,l,param[[3]][[l]],N,delta,grid)[,,j]		
						sig.comp <- f.comp %*% solve(g.comp) %*% h.comp  
					}
					frf.sig <- frf.sig + sig.comp
				}
				delta <- sigex.delta(mdl,noisecomps)
				if(Mod(sum(delta*exp(-1i*seq(0,length(delta)-1)*(j-1)*pi/grid))) < 10^(-8))
				{
					delta <- sigex.delta(mdl,sigcomps)
					frf.wk[,,j] <- (frf.sig %*% f.noise[,,j] + (diag(N) - frf.sig) %*% f.sig[,,j])*
						Mod(sum(delta*exp(-1i*seq(0,length(delta)-1)*(j-1)*pi/grid)))^{-2}
				}
				delta <- sigex.delta(mdl,sigcomps)
				if(Mod(sum(delta*exp(-1i*seq(0,length(delta)-1)*(j-1)*pi/grid))) < 10^(-8))
				{
					delta <- sigex.delta(mdl,noisecomps)
					frf.wk[,,j] <- (frf.sig %*% f.noise[,,j] + (diag(N) - frf.sig) %*% f.sig[,,j])*
						Mod(sum(delta*exp(-1i*seq(0,length(delta)-1)*(j-1)*pi/grid)))^{-2}
				}
			}  
		}
	 	if(!flag.zero) { frf.wk[,,j] <- f.sig[,,j] %*% solve(f.all[,,j]) %*% f.noise[,,j] }
	}	

	return(frf.wk)
}
 