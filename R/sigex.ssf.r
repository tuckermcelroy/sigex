sigex.ssf <- function(data,param,mdl,sigcomps)
{

	###############################
	#   sigex.ssf
	#	by Tucker McElroy
	#
	#	Puts model into SSF commpatible with Assumption A.
	#		runs ssf filter and smoother, properly initialized.
	#		param must be in format yielded by sigex.default
	#		sigcomps provides indices of the desired components
	#     Output is the signal estimates (TxN matrix) and precisions (TxN matrix)
	#		as a list object.
	#	NOTE: this version presumes diff'd components are white noise
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# put mdl in SSF
	agg.vec <- NULL
	phi.mat <- NULL
	for(i in 1:length(mdl[[3]]))
	{
		delta.comp <- mdl[[3]][[i]]
		len.delta <- length(delta.comp)-1
		if(len.delta == 0)
		{
			agg.vec <- cbind(agg.vec,diag(N))
			if(i==1) { phi.mat <- matrix(0,N,N) } else {
			phi.mat <- cbind(phi.mat,matrix(0,nrow=dim(phi.mat)[1],ncol=N))
			phi.mat <- rbind(phi.mat,matrix(0,nrow=N,ncol=dim(phi.mat)[2]))
			}
		}
		if(len.delta == 1)
		{
			agg.vec <- cbind(agg.vec,diag(N))
			new.phi <- -1*(delta.comp[-1]/delta.comp[1])*diag(N)
			if(i==1) { phi.mat <- new.phi } else {
			phi.row <- dim(phi.mat)[1]
			phi.col <- dim(phi.mat)[2]
			phi.mat <- cbind(phi.mat,matrix(0,nrow=phi.row,ncol=N))
			phi.mat <- rbind(phi.mat,matrix(0,nrow=N,ncol=phi.col),new.phi)
			}
		}
		if(len.delta > 1)
		{
			agg.vec <- cbind(agg.vec,diag(N),matrix(0,N,N*(len.delta-1)))
			new.phi <- -1*delta.comp[-1]/delta.comp[1]
			new.phi <- rbind(new.phi,diag(len.delta)[1:(len.delta-1),])
			new.phi <- new.phi %x% diag(N)
			if(i==1) { phi.mat <- new.phi } else {
			phi.row <- dim(phi.mat)[1]
			phi.col <- dim(phi.mat)[2]
			phi.mat <- cbind(phi.mat,matrix(0,nrow=phi.row,ncol=dim(new.phi)[2]))
			phi.mat <- rbind(phi.mat,cbind(matrix(0,nrow=dim(new.phi)[1],ncol=phi.col),new.phi))
			}
		}
	}

#	print(phi.mat)
#	print(agg.vec)

	# initialize Kalman filter for data up to time d+2, where d is order of delta(B)
	d.delta <- length(sigex.delta(mdl,0))-1
	T.init <- d.delta+2
	mdl.short <- mdl
	for(k in 1:N) { mdl.short[[4]][[k]] <- mdl[[4]][[k]][1:T.init,] }
	init.data <- data[1:T.init,]
	state.init <- NULL
	for(i in 1:length(mdl[[3]]))
	{
		delta.comp <- mdl[[3]][[i]]
		len.delta <- length(delta.comp)-1
		filters.array <- sigex.signal(init.data,param,mdl,i)
		if(i==sigex.whichtrend(mdl)) { meanFlag <- 1 } else { meanFlag <- 0 }
		signals.array <- sigex.extract(init.data,filters.array,mdl.short,param,meanFlag)
		if(len.delta == 0) {
			state.init <- rbind(state.init,signals.array[[1]][T.init,])
		} else {
			state.init <- rbind(state.init,signals.array[[1]][T.init:(T.init-len.delta+1),])
		}
	}
	state.init <- matrix(t(state.init),ncol=1)
#	var.init

	# given initial state estimate and precision, run state space filter
	state.len <- length(state.init)
	state.vecs <- state.init
	state.vars <- array(0,c(state.len,state.len,(T-T.init+1)))
	state.vars[,,1] <- var.init
	pred.vecs <- NULL
	pred.vars <- array(0,c(state.len,state.len,(T-T.init)))
	resid.vecs <- NULL
	resid.vars <- array(0,c(N,N,(T-T.init)))
	for(t in (T.init+1):T)
	{
		pred.vec <- phi.mat %*% state.vecs[,(t-T.init)]
		pred.var <- phi.mat %*% state.vars[,,(t-T.init)] %*% t(phi.mat)
		resid.vec <- data[t,] - agg.vec %*% pred.vec
		resid.var <- agg.vec %*% pred.var %*% t(agg.vec)
		kalman <- pred.var %*% agg.vec %*% solve(resid.var)
		state.vec <- pred.vec + kalman %*% resid.vec
		state.var <- (diag(length(state.init)) - kalman %*% t(agg.mat)) %*% pred.var
		pred.vecs <- cbind(preds.vec,pred.vec)
		pred.vars[,,(t-T.init)] <- pred.var
		resid.vecs <- cbind(resid.vecs,resid.vec)
		resid.vars[,,(t-T.init)] <- resid.var
		state.vecs <- cbind(state.vecs,state.vec)
		state.vars[,,(t-T.init+1)] <- state.var
	}

	# run state space smoother
	# HERE
	smooth.vecs <- state.vecs[,(T-T.init+1)]
	smooth.vars <- array(0,c(state.len,state.len,T))
	smooth.vars[,,T] <- state.vars[,,(T-T.init+1)]
	for(t in T:2)
	{
		# JL - ERROR on line 122 just commented it out
	  #smooth.vec <- 	Jmat %*% (smooth.vecs[,(T-t+1)] - pred.vec[,

		smooth.vecs <- cbind(smooth.vec,smooth.vecs)
	}

}
















