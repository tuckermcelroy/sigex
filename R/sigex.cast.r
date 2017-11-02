sigex.cast <- function(psi,mdl,data,leads,plotit=TRUE)
{

	###########################
	#	sigex.cast
	#		by Tucker McElroy
	#
	#   Computes forecasts/aftcasts of a latent multivariate component process,
	#		where regression effects have been removed
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#	mdl is the model structure
	#	x is the data, NxT matrix
	#	leads is an integer sequence of desired forecasts/aftcasts; include index i
	#		to obtain an estimate of x_i.  So 1<=i<=T corresponds to data,
	#		and i<1 to aftcasts, i>T to forecasts
	#	returns Nxlength(leads) matrix of casts; optional plotting
	#		
	##############################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	psi <- Re(psi)

	indices <- union(seq(1,T),leads)
	aft.index <- min(indices)
	fore.index <- max(indices)

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	acf.mat <- matrix(0,nrow=N*length(indices),ncol=N)
	
	# get xi portion
	ind <- 0
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	for(i in 1:length(mdl[[3]]))
	{
		vrank <- mdl[[1]][[i]]
		D.dim <- length(vrank)
		L.dim <- sum(A.mat[,as.vector(vrank)])
		L.psi <- NULL
		if(L.dim > 0) L.psi <- psi[(ind+1):(ind+L.dim)]
		ind <- ind+L.dim
		D.psi <- psi[(ind+1):(ind+D.dim)]
		ind <- ind+D.dim
		L.mat <- sigex.param2gcd(L.psi,N,as.vector(vrank))
		L.par[[i]] <- L.mat
		D.par[[i]] <- D.psi
	}

	# get beta portion
	beta.len <- 0
	for(i in 1:N) 
	{
		beta.len <- beta.len + dim(mdl[[4]][[i]])[2]
	}
	beta.par <- as.vector(psi[(length(psi)-beta.len+1):length(psi)])

	# get zeta portion
	if(length(psi)-beta.len-ind > 0) {
		zeta <- psi[(ind+1):(length(psi)-beta.len)] }
	ind <- 0
	for(i in 1:length(mdl[[3]]))
	{
		mdlType <- mdl[[2]][i]	
		delta <- mdl[[3]][[i]]
		zetalen <- sigex.zetalen(mdlType)
		if(zetalen > 0) {
			subzeta <- zeta[(ind+1):(ind+zetalen)]
			zeta.par[[i]] <- sigex.zeta2par(subzeta,mdlType,delta,N)
		}
		ind <- ind + zetalen
	
		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			zeta.par[[i]],N,delta,length(indices))		
	}

	x.acf <- array(acf.mat,dim=c(N,length(indices),N))
	reg.vec <- beta.par	

	# subtract regression effects
	ind <- 0
	data.diff <- data
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		data.diff[,k] <- data.diff[,k] - reg.mat %*% reg.vec[(ind+1):(ind+len)]
		ind <- ind+len
	}

	# difference the data
	fulldiff <-  sigex.delta(mdl,0)
	del <- length(fulldiff) - 1
	x.diff <- as.matrix(filter(data.diff,fulldiff,method="convolution",
		sides=1)[length(fulldiff):T,])
	Tdiff <- dim(x.diff)[1]
	x.diff <- t(x.diff)
 
	fore.cast <- NULL
	aft.cast <- NULL
	if(fore.index > T)
	{
		x.fore <- cbind(x.diff,matrix(1i,N,(fore.index-T)))
#		diff.cast <- mvar.cast(x.acf,x.fore)
		diff.cast <- mvar.forecast(x.acf,x.fore)[[1]]
		## add forecasted regressors???
		if(del > 0) {
			fore.cast <- as.matrix(filter(init = matrix(data.diff[del:1,],ncol=N),
				x=t(diff.cast)/fulldiff[1],filter=-1*fulldiff[-1]/fulldiff[1],
				method="recursive"))
		} else { fore.cast <- t(diff.cast) }
		fore.cast <- as.matrix(fore.cast[(Tdiff+1):(Tdiff+fore.index-T),])
		fore.cast <- t(fore.cast)		
	}
	if(aft.index < 1)
	{
		x.rev <- t(as.matrix(t(x.diff)[seq(Tdiff,1),]))
		x.aft <- cbind(x.rev,matrix(1i,N,(1-aft.index)))
#		diff.cast <- mvar.cast(aperm(x.acf,c(3,2,1)),x.aft) 
		diff.cast <- mvar.forecast(aperm(x.acf,c(3,2,1)),x.aft)[[1]]
 		## add forecasted regressors???
		if(del > 0) {
			aft.cast <- as.matrix(filter(init = matrix(data.diff[(Tdiff+1):T,],ncol=N),
				x=t(diff.cast)/fulldiff[del+1],filter=-1*rev(fulldiff)[-1]/fulldiff[del+1],
				method="recursive"))
		} else { aft.cast <- t(diff.cast) }
		aft.cast <- as.matrix(aft.cast[(Tdiff+1):(Tdiff+1-aft.index),])
		aft.cast <- as.matrix(aft.cast[seq(1-aft.index,1),])
		aft.cast <- t(aft.cast)		
	}
	x.casted <- cbind(aft.cast,t(data.diff),fore.cast)
	x.real <- x.casted
	if(length(aft.cast) > 0) { x.real[,1:dim(aft.cast)[2]] <- NA }
	if(length(fore.cast) > 0) { x.real[,(dim(x.real)[2]-dim(fore.cast)[2]):dim(x.real)[2]] <- NA }

	if(plotit && (N <= 6)) {
		if(N <= 4) { par(mfrow=c(N,1)) } 
		if((5 <= N) && (N <= 6)) { par(mfrow=c(3,2)) } 
	      start.date <- day2date(start(data)[2]-2+aft.index,c(1,1,start(data)[1]))
		start.day <- date2day(start.date[1],start.date[2],start.date[3])
		begin <- c(start.date[3],start.day)
		for(i in 1:N) 
		{
			plot(ts(x.casted[i,],start=begin,frequency=frequency(data)),ylab="",xlab="Year")
			lines(ts(x.real[i,],start=begin,frequency=frequency(data)),col=2)
		}
	 }

	return(x.casted)
}


