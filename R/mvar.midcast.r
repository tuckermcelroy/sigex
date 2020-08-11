mvar.midcast <- function(x.acf,z,delta,debug=FALSE)
{

	##########################################################################
	#
	#	mvar.midcast
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
	#	Purpose: compute multi-step imputations and predictors of a multivariate process
	#		via Levinson-Durbin algorithm with  missing values
	#	Background:
	#		A multivariate difference-stationary process x_t with
	#			w_t = delta(B) x_t
	#		may be observed with missing values, and one wants to compute
	#		Gaussian conditional expectations of missing values (midcasts),
	#		or future values (forecasts), or past values (aftcasts).
	#		Also of interest is the Gaussian likelihood resulting from
	#		such a sample, and the residuals.
	#		It is required that at least d
	#		contiguous values be observed, where d is the order of delta(B).
	#		If the first d values are contiguous, we can do a forward pass;
	#		otherwise, the first set of d contiguous values starts after
	#		time index 1, and is marked by t.hash.  In this case, we must
	#		also do a backward pass, involving aftcasts.
	#	Inputs:
	#		x.acf: array of dimension N x T x N of autocovariances for process w_t,
	#			where there are N series, of total length T each.
	#		z: raw data as N x T matrix with missing values at various time points.
	#			Missing values are at any 1 <= t <= T, and occur for some of the N series
  #     (the ragged case), and are denoted with a 1i.  That is,
	#			Im(z[,t]) = rep(1i,N) or subset thereof encodes missing values.
	#		delta: differencing polynomial (corresponds to delta(B) in Background)
	#			written in format c(delta0,delta1,...,deltad)
  #   debug: set to TRUE if lik values should be printed to screen
	#	Notes: to get H forecasts, append matrix(1i,N,H) to input x.  To get aftcasts,
	#		prepend the same.  T will be the second dimension of z, and includes
	#		the spots taken by aftcasts and forecasts.  (So the T for the original
	#		dataset could be less than the T used in this function.)
 	#	Outputs:
	#		list containing casts.x, casts.var, c(Qseq,logdet), and eps
	#		casts.x: N x H matrix of backcasts, midcasts, aftcasts, where H
	#			is the total number of time indices with missing values.
	#			So if times.na is a subset of seq(1,T) corresponding to indices
	#			with missing values, we can fill in all NAs via	z[,times.na] <- casts.x.
  #     If a subset is missing at some time t, casts.x for that time t contains
  #     the casts together with the known values.
	#		casts.var: NH x NH matrix of covariances of casting errors.
	#			note that casts.var.array <- array(casts.var,c(N,H,N,H))
	#			corresponds to cast.var.array[,j,,k] equal to the
	#			covariance between the jth and kth casting errors.
  #     If a subset is missing at some time t, the corresponding block of casts.var
  #     will have zeros corresponding to the known values.
	#		Qseq: quadratic form portion of the Gaussian divergence based on
	#			missing value formulation (McElroy and Monsell, 2015 JASA)
	#		logdet: log determinant portion of the Gaussian divergence
	#		eps: residuals from casting recursions, defined in the manner
	#			discussed in Casting paper.  Dimension is N x T-(H+d)
	#
	####################################################################

	thresh <- 10^(-16)
	N <- dim(z)[1]
	T <- dim(z)[2]
	all.series <- seq(1,N)
	all.indices <- seq(1,T)
	full.indices <- all.indices[colSums(Im(z)==1)==0]
	cast.indices <- setdiff(all.indices,full.indices)
	ragged <- list()
	leads.rag <- NULL
	for(t in 1:length(cast.indices))
	{
	  rag.series <- all.series[Im(z[,cast.indices[t]])==1,drop=FALSE]
	  if(length(rag.series)<=N)
	  {
	    ragged[[length(ragged)+1]] <- rag.series
      leads.rag <- c(leads.rag,cast.indices[t])
	  }
	}
	d <- length(delta) - 1

	# This version does not presume that the first d values are not missing
	# Find t.hash, earliest time for which d contiguous values follow,
	#	such that data at times t.hash+1,...,t.hash+d exists
	t.hash <- NULL
	ind.data <- c(0,cast.indices,(T+1))
	gaps <- c(diff(ind.data),0)
	if(max(gaps) > d) { t.hash <- min(ind.data[gaps > d]) }

	if(d > 0) {
	# t=t.hash+d case as initialization
	#  get predictors based on observations t.hash+1:t
	l.pred <- t(-1*delta[1]^{-1}*delta[(d+1):2] %x% diag(N))
	l.derp <- t(-1*delta[d+1]^{-1}*delta[d:1] %x% diag(N))
	v.pred <- as.matrix(delta[1]^{-2}*x.acf[,1,])
	v.derp <- as.matrix(delta[d+1]^{-2}*x.acf[,1,])
	#  get casts and covars of observations t.hash+1:t based on sigma-field_{t.hash+1:t}
	#   Note: store more than needed in preds.x, makes it easier for indexing later
	preds.x <- Re(z[,1:(t.hash+d),drop=FALSE])
	new.covar <- NULL
	casts.x <- NULL
	casts.var <- NULL
	eps <- NULL
	Qseq <- 0
	logdet <- 0
	# track indices of casted variables up to present time
	cast.index.t <- intersect(cast.indices,seq(t.hash+1,t.hash+d))
	t.star <- T
	t.len <- d

	# Forward Pass:
	if(t.hash < T-d) {
	for(t in (t.hash+d+1):T)
	{
	  # determine whether full info, or partial/completely missing
	  #  base case: full info
	  select.mat <- diag(N)
	  omit.mat <- NULL
	  raggeds <- NULL
	  rags.ind <- leads.rag %in% t
	  if(sum(rags.ind)>0) # partial/completely missing
	  {
	    raggeds <- ragged[[seq(1,length(leads.rag))[rags.ind]]]
	    if(length(raggeds)<N) # partial missing
	    {
  	    select.mat <- diag(N)[-raggeds,,drop=FALSE]
	      omit.mat <- diag(N)[raggeds,,drop=FALSE]
	    } else # completely missing
	    {
	      select.mat <- NULL
	      omit.mat <- diag(N)
	    }
	  }
	  non.raggeds <- setdiff(seq(1,N),raggeds)

		#  get casts and covars of observations t.hash+1:t based on sigma-field_{t.hash+1:t}
		# first, construct preds.x from known (unstored) entries and stored casts;
	  #   preds.x is E [ X_{t-1} | F_{t-1} ]
		if(length(cast.index.t) > 0) { preds.x[,cast.index.t] <- casts.x }
	  # obtain E [ x_t | F_{t-1} ]
		new.pred <- l.pred %*% matrix(preds.x[,(t-t.len):(t-1)],ncol=1)

	  # second, get prediction variance, to obtain new.var given by Var [ x_t | F_{t-1} ]
		#   cast.index.tlen tracks cast indices up to now, looking back t.len time points,
		#   where t.len is the range at which all dependence is effectively nil
		cast.index.tlen <- intersect(cast.indices,seq(t-t.len,t-1))
		cast.len <- length(cast.index.tlen)
		if(cast.len==0) # no casts within t.len time points
		{
		  new.var <- v.pred
		} else # at least one cast within t.len time points
		{
		  casts.var.array <- array(casts.var,c(N,length(cast.index.t),N,length(cast.index.t)))
		  range.t <- (length(cast.index.t)-cast.len+1):length(cast.index.t)
		  casts.var.array <- casts.var.array[,,,range.t,drop=FALSE]
		  if(t-1-t.len>0) { l.pred <- cbind(matrix(0,N,N*(t-1-t.len)),l.pred) }
		  l.array <- array(l.pred,c(N,N,t-1))
		  l.array <- l.array[,,cast.index.tlen,drop=FALSE]
		  l.pred.tlen <- matrix(l.array,nrow=N)
			new.var <- v.pred + l.pred.tlen %*%
		    matrix(casts.var.array[,range.t,,,drop=FALSE],nrow=cast.len*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		}

		# third, update casts.x by changing the stored portions of
		#   E [ X_{t-1} | F_{t-1} ] to E [ X_{t-1} | F_t ] and
		#   appending E [ x_t | F_t ] if partially/completely missing
    if(cast.len>0)  # at least one cast within t.len time points
    {
      if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
      {
		    new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		    update <- matrix(new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred),nrow=N)
		    casts.x <- casts.x + update
      }
    }
    if(length(raggeds)>0)   # add new cast E [ x_t | F_t ] if partially/completely missing
    {
      new.cast <- Re(z[,t,drop=FALSE])
      partial.cast <- omit.mat %*% new.pred
      if(length(raggeds)<N)  # partial missing case
      {
        partial.cast <- partial.cast + omit.mat %*% new.var %*% t(select.mat) %*%
          solve(select.mat %*% new.var %*% t(select.mat)) %*%
          (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
      }
      new.cast[raggeds,1] <- partial.cast
      casts.x <- cbind(casts.x,new.cast)
    }

    # fourth, update casts.var by changing the stored portions of
		#   Var [ X_{t-1} | F_{t-1} ] to Var [ X_{t-1} | F_t ] and
		#   appending new covariances if partially/completely missing
		if(cast.len>0)  # at least one cast within t.len time points
		{
		  new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		  if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
		  {
		    update <- new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      select.mat %*% t(new.covar)
		    casts.var <- casts.var - update
		  }
		}
		if(length(raggeds)>0)   # add new covar [ x_t | F_t ] if partially/completely missing
		{
		  proj <- diag(N) # completely missing case
		  if(length(raggeds)<N)  # partial missing case
		  {
		    proj <- diag(N) - t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      select.mat %*% new.var
		  }
		  # now augment casts.var
		  # special case: no cast within t.len time points.
		  #   Either: (i) this is the first cast,
		  #   Or: (ii) previous casts were long ago
		  if(cast.len==0) # special case
		  {
		    if(length(casts.var)==0)  # (i) of special case
		    {
		      new.block <- NULL
		    } else  # (ii) of special case
		    {
		      new.covar <- matrix(0,nrow=length(cast.index.t)*N,ncol=N)
		      new.block <- new.covar %*% proj
		    }
		  } else
		  {
		    new.block <- new.covar %*% proj
		  }
		  # do regular case (and special case) augmentation
		  casts.var <- rbind(cbind(casts.var,new.block),
		                     t(rbind(new.block,t(new.var %*% proj))))
		}

	  # fifth, get ragged residuals
		if(length(raggeds)>0)  # case of partially/completely missing
		{
		  new.eps <- matrix(rep(1i,N),ncol=1)
		  new.det <- 1
		  if(length(raggeds)<N)  # partial missing case
		  {
		    partial.eps <- solve(t(chol(select.mat %*% new.var %*% t(select.mat)))) %*%
		      (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
		    new.eps[non.raggeds,1] <- partial.eps
		    new.det <- det(chol(select.mat %*% new.var %*% t(select.mat)))
		  }
		} else # case of full info
		{
		  new.eps <- solve(t(chol(new.var))) %*% (Re(z[,t,drop=FALSE]) - new.pred)
		  new.det <- det(new.var)
		}
		eps <- rbind(eps,new.eps)
		Qseq <- Qseq + t(Re(new.eps)) %*% Re(new.eps)
		logdet <- logdet + log(new.det)

		# updating
		cast.index.t <- intersect(cast.indices,seq(t.hash+1,t))
		preds.x <- Re(z[,1:t,drop=FALSE])

		#  get fore and aft predictors based on observations (t.hash+1):t
		#  Notes: in non-stationary case, we obtain predictors for time t+1,
		#   based on data over times t.hash+1,...,t.hash+d,...,t.
		#   That is, the predictors computed below are used in the next
		#   increment of the t-loop.
		if(t==(t.hash+d+1))
		{

			u.seq <- solve(x.acf[,1,]) %*% x.acf[,2,]
			l.seq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			gam.Seq <- x.acf[,2,]
			gam.Flip <- x.acf[,2,]
			c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
			d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
			a.seq <- (delta[(d+1):1] %x% diag(N)) %*% u.seq
			b.seq <- (delta[(d+1):1] %x% diag(N)) %*% l.seq

			u.qes <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			l.qes <- solve(x.acf[,1,]) %*% x.acf[,2,]
			gam.qeS <- t(x.acf[,2,])
			gam.pilF <- t(x.acf[,2,])
			c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
			d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
			a.qes <- (delta[(d+1):1] %x% diag(N)) %*% u.qes
			b.qes <- (delta[(d+1):1] %x% diag(N)) %*% l.qes

		} else # case of t.hash+d+2 <= t <= T
		{
			if(t.star == T)
			{

				pacf <- x.acf[,t+1-d-t.hash,] - gam.Seq %*% u.seq
				l.factor <- solve(c.mat) %*% t(pacf)
				new.l <- l.seq - u.seq %*% l.factor
				u.factor <- solve(d.mat) %*% pacf
				new.u <- u.seq - l.seq %*% u.factor
				l.seq <- rbind(l.factor,new.l)
				u.seq <- rbind(new.u,u.factor)
				gam.Seq <- cbind(x.acf[,t+1-d-t.hash,],gam.Seq)
				gam.Flip <- rbind(gam.Flip,x.acf[,t+1-d-t.hash,])
				c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
				d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
				a.next <- a.seq - b.seq %*% u.factor
				b.next <- b.seq - a.seq %*% l.factor
				a.seq <- rbind(a.next,0*diag(N)) + (c(rep(0,t-d-1-t.hash),delta[(d+1):1]) %x% diag(N)) %*% u.factor
				b.seq <- rbind(0*diag(N),b.next) + (c(delta[(d+1):1],rep(0,t-d-1-t.hash)) %x% diag(N)) %*% l.factor

				fcap <- t(x.acf[,t+1-d-t.hash,]) - gam.qeS %*% u.qes
				l.rotcaf <- solve(c.tam) %*% t(fcap)
				wen.l <- l.qes - u.qes %*% l.rotcaf
				u.rotcaf <- solve(d.tam) %*% fcap
				wen.u <- u.qes - l.qes %*% u.rotcaf
				l.qes <- rbind(l.rotcaf,wen.l)
				u.qes <- rbind(wen.u,u.rotcaf)
				gam.qeS <- cbind(t(x.acf[,t+1-d-t.hash,]),gam.qeS)
				gam.pilF <- rbind(gam.pilF,t(x.acf[,t+1-d-t.hash,]))
				c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
				d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
				a.next <- a.qes - b.qes %*% u.rotcaf
				b.next <- b.qes - a.qes %*% l.rotcaf
				a.qes <- rbind(0*diag(N),a.next) + (c(delta[(d+1):1],rep(0,t-d-1-t.hash)) %x% diag(N)) %*% u.rotcaf
				b.qes <- rbind(b.next,0*diag(N)) + (c(rep(0,t-d-1-t.hash),delta[(d+1):1]) %x% diag(N)) %*% l.rotcaf

			}
			if(sqrt(sum(diag(l.factor %*% t(l.factor)))) < thresh)
			{
				t.star <- min(t.star,t)
			}
		}
		t.len <- dim(b.seq)[1]/N
		l.pred <- t(-1*delta[1]^{-1}*c(rep(0,t.len-d),delta[(d+1):2]) %x% diag(N))
		l.derp <- t(-1*delta[d+1]^{-1}*c(delta[d:1],rep(0,t.len-d)) %x% diag(N))
		l.pred <- l.pred + delta[1]^{-1}*t(b.seq)
		l.derp <- l.derp + delta[d+1]^{-1}*t(b.qes)
		v.pred <- delta[1]^{-2}*(x.acf[,1,] - matrix(aperm(x.acf[,(t.len-d+1):2,,drop=FALSE],c(1,3,2)),nrow=N) %*% l.seq)
		v.derp <- delta[d+1]^{-2}*(x.acf[,1,] - matrix(aperm(x.acf[,(t.len-d+1):2,,drop=FALSE],c(3,1,2)),nrow=N) %*% l.qes)
	} }

	# Backward Pass:
	if(t.hash > 0) {
	for(t in (t.hash):1)
	{
	  # determine whether full info, or partial/completely missing
	  #  base case: full info
	  select.mat <- diag(N)
	  omit.mat <- NULL
	  raggeds <- NULL
	  rags.ind <- leads.rag %in% t
	  if(sum(rags.ind)>0) # partial/completely missing
	  {
	    raggeds <- ragged[[seq(1,length(leads.rag))[rags.ind]]]
	    if(length(raggeds)<N) # partial missing
	    {
	      select.mat <- diag(N)[-raggeds,,drop=FALSE]
	      omit.mat <- diag(N)[raggeds,,drop=FALSE]
	    } else # completely missing
	    {
	      select.mat <- NULL
	      omit.mat <- diag(N)
	    }
	  }
	  non.raggeds <- setdiff(seq(1,N),raggeds)

	  #  get casts and covars of observations t:T based on sigma-field_{t:T}
	  # first, construct preds.x from known (unstored) entries and stored casts;
	  #   preds.x is E [ X_{t-1} | F_{t+1} ]
	  if(length(cast.index.t) > 0) { preds.x[,cast.index.t] <- casts.x }
	  # obtain E [ x_t | F_{t-1} ]
		new.pred <- l.derp %*% matrix(preds.x[,(t+1):(t+t.len)],ncol=1)

		# second, get prediction variance, to obtain new.var given by Var [ x_t | F_{t+1} ]
		#   cast.index.tlen tracks cast indices up to now, looking forward t.len time points,
		#   where t.len is the range at which all dependence is effectively nil
		cast.index.tlen <- intersect(cast.indices,seq(t+1,t+t.len))
		cast.len <- length(cast.index.tlen)
		if(cast.len==0) # no casts within t.len time points
		{
		  new.var <- v.derp
		} else # at least one cast within t.len time points
		{
		  casts.var.array <- array(casts.var,c(N,length(cast.index.t),N,length(cast.index.t)))
#		  range.t <- (length(cast.index.t)-cast.len+1):length(cast.index.t)
		  range.t <- 1:cast.len
		  casts.var.array <- casts.var.array[,,,range.t,drop=FALSE]
		  l.array <- array(l.derp,c(N,N,T-t))
		  l.array <- l.array[,,cast.index.tlen-t,drop=FALSE]
		  l.pred.tlen <- matrix(l.array,nrow=N)
		  new.var <- v.derp + l.pred.tlen %*%
		    matrix(casts.var.array[,range.t,,,drop=FALSE],nrow=cast.len*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		}

		# third, update casts.x by changing the stored portions of
		#   E [ X_{t+1} | F_{t+1} ] to E [ X_{t+1} | F_t ] and
		#   appending E [ x_t | F_t ] if partially/completely missing
		if(cast.len>0)  # at least one cast within t.len time points
		{
		  if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
		  {
		    new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		    update <- matrix(new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		                       (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred),nrow=N)
		    casts.x <- casts.x + update
		  }
		}
		if(length(raggeds)>0)   # add new cast E [ x_t | F_t ] if partially/completely missing
		{
		  new.cast <- Re(z[,t,drop=FALSE])
		  partial.cast <- omit.mat %*% new.pred
		  if(length(raggeds)<N)  # partial missing case
		  {
		    partial.cast <- partial.cast + omit.mat %*% new.var %*% t(select.mat) %*%
		      solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
		  }
		  new.cast[raggeds,1] <- partial.cast
		  casts.x <- cbind(new.cast,casts.x)
		}

		# fourth, update casts.var by changing the stored portions of
		#   Var [ X_{t+1} | F_{t+1} ] to Var [ X_{t+1} | F_t ] and
		#   appending new covariances if partially/completely missing
		if(cast.len>0)  # at least one cast within t.len time points
		{
		  new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
		  if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
		  {
		    update <- new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      select.mat %*% t(new.covar)
		    casts.var <- casts.var - update
		  }
		}
		if(length(raggeds)>0)   # add new covar [ x_t | F_t ] if partially/completely missing
		{
		  proj <- diag(N) # completely missing case
		  if(length(raggeds)<N)  # partial missing case
		  {
		    proj <- diag(N) - t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
		      select.mat %*% new.var
		  }
		  # now augment casts.var
		  # special case: no cast within t.len time points.
		  #   Either: (i) this is the first cast,
		  #   Or: (ii) previous casts were long ago
		  if(cast.len==0) # special case
		  {
		    if(length(casts.var)==0)  # (i) of special case
		    {
		      new.block <- NULL
		    } else  # (ii) of special case
		    {
		      new.covar <- matrix(0,nrow=length(cast.index.t)*N,ncol=N)
		      new.block <- new.covar %*% proj
		    }
		  } else
		  {
		    new.block <- new.covar %*% proj
		  }
		  # do regular case (and special case) augmentation
		  casts.var <- rbind(t(rbind(t(new.var %*% proj),new.block)),
	  	                   cbind(new.block,casts.var))
		}

		# fifth, get ragged residuals
		if(length(raggeds)>0)  # case of partially/completely missing
		{
		  new.eps <- matrix(rep(1i,N),ncol=1)
		  new.det <- 1
		  if(length(raggeds)<N)  # partial missing case
		  {
		    partial.eps <- solve(t(chol(select.mat %*% new.var %*% t(select.mat)))) %*%
		      (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
		    new.eps[non.raggeds,1] <- partial.eps
		    new.det <- det(chol(select.mat %*% new.var %*% t(select.mat)))
		  }
		} else # case of full info
		{
		  new.eps <- solve(t(chol(new.var))) %*% (Re(z[,t,drop=FALSE]) - new.pred)
		  new.det <- det(new.var)
		}
		eps <- rbind(eps,new.eps)
		Qseq <- Qseq + t(Re(new.eps)) %*% Re(new.eps)
		logdet <- logdet + log(new.det)

		# updating
		cast.index.t <- intersect(cast.indices,seq(t,T))
		preds.x <- Re(z[,1:T,drop=FALSE])

		#  get predictors based on observations t+1:T
		#  Notes: in non-stationary case, we obtain predictors for time t-1,
		#   based on data over times t,...,T.
		#   That is, the predictors computed below are used in the next
		#   increment of the t-loop.
		if(t==(T-d))	# if no forward pass happened, initialize
		{

			u.seq <- solve(x.acf[,1,]) %*% x.acf[,2,]
			l.seq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			gam.Seq <- x.acf[,2,]
			gam.Flip <- x.acf[,2,]
			c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
			d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
			a.seq <- (delta[(d+1):1] %x% diag(N)) %*% u.seq
			b.seq <- (delta[(d+1):1] %x% diag(N)) %*% l.seq

			u.qes <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			l.qes <- solve(x.acf[,1,]) %*% x.acf[,2,]
			gam.qeS <- t(x.acf[,2,])
			gam.pilF <- t(x.acf[,2,])
			c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
			d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
			a.qes <- (delta[(d+1):1] %x% diag(N)) %*% u.qes
			b.qes <- (delta[(d+1):1] %x% diag(N)) %*% l.qes

		} else # case of 1 <= t <= t.hash and t.hash < T-d
		{
			if(t.star == T)
			{

				pacf <- x.acf[,T-d-t+2,] - gam.Seq %*% u.seq
				l.factor <- solve(c.mat) %*% t(pacf)
				new.l <- l.seq - u.seq %*% l.factor
				u.factor <- solve(d.mat) %*% pacf
				new.u <- u.seq - l.seq %*% u.factor
				l.seq <- rbind(l.factor,new.l)
				u.seq <- rbind(new.u,u.factor)
				gam.Seq <- cbind(x.acf[,T-d-t+2,],gam.Seq)
				gam.Flip <- rbind(gam.Flip,x.acf[,T-d-t+2,])
				c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
				d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
				a.next <- a.seq - b.seq %*% u.factor
				b.next <- b.seq - a.seq %*% l.factor
				a.seq <- rbind(a.next,0*diag(N)) + (c(rep(0,T-d-t),delta[(d+1):1]) %x% diag(N)) %*% u.factor
				b.seq <- rbind(0*diag(N),b.next) + (c(delta[(d+1):1],rep(0,T-d-t)) %x% diag(N)) %*% l.factor

				fcap <- t(x.acf[,T-d-t+2,]) - gam.qeS %*% u.qes
				l.rotcaf <- solve(c.tam) %*% t(fcap)
				wen.l <- l.qes - u.qes %*% l.rotcaf
				u.rotcaf <- solve(d.tam) %*% fcap
				wen.u <- u.qes - l.qes %*% u.rotcaf
				l.qes <- rbind(l.rotcaf,wen.l)
				u.qes <- rbind(wen.u,u.rotcaf)
				gam.qeS <- cbind(t(x.acf[,T-d-t+2,]),gam.qeS)
				gam.pilF <- rbind(gam.pilF,t(x.acf[,T-d-t+2,]))
				c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
				d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
				a.next <- a.qes - b.qes %*% u.rotcaf
				b.next <- b.qes - a.qes %*% l.rotcaf
				a.qes <- rbind(0*diag(N),a.next) + (c(delta[(d+1):1],rep(0,T-d-t)) %x% diag(N)) %*% u.rotcaf
				b.qes <- rbind(b.next,0*diag(N)) + (c(rep(0,T-d-t),delta[(d+1):1]) %x% diag(N)) %*% l.rotcaf

			}
			if(sqrt(sum(diag(l.factor %*% t(l.factor)))) < thresh)
			{
				t.star <- min(t.star,t)
			}
		}
		t.len <- dim(b.seq)[1]/N
		l.pred <- t(-1*delta[1]^{-1}*c(rep(0,t.len-d),delta[(d+1):2]) %x% diag(N))
		l.derp <- t(-1*delta[d+1]^{-1}*c(delta[d:1],rep(0,t.len-d)) %x% diag(N))
		l.pred <- l.pred + delta[1]^{-1}*t(b.seq)
		l.derp <- l.derp + delta[d+1]^{-1}*t(b.qes)
		v.pred <- delta[1]^{-2}*(x.acf[,1,] - matrix(aperm(x.acf[,(t.len-d+1):2,,drop=FALSE],c(1,3,2)),nrow=N) %*% l.seq)
		v.derp <- delta[d+1]^{-2}*(x.acf[,1,] - matrix(aperm(x.acf[,(t.len-d+1):2,,drop=FALSE],c(3,1,2)),nrow=N) %*% l.qes)
	} }

	} else  # d = 0
	{

	# initializations based on no data
	new.pred <- matrix(0,nrow=N,ncol=1)
	v.pred <- as.matrix(x.acf[,1,])
	new.covar <- NULL
	casts.x <- NULL
	casts.var <- NULL
	eps <- NULL
	Qseq <- 0
	logdet <- 0
	cast.index.t <- NULL
	cast.index.tlen <- NULL
	t.star <- T
	t.len <- 0

	# Forward Pass:
	if(t.hash < T) {
  for(t in (t.hash+1):T)
	{
	  # determine whether full info, or partial/completely missing
	  #  base case: full info
	  select.mat <- diag(N)
	  omit.mat <- NULL
	  raggeds <- NULL
	  rags.ind <- leads.rag %in% t
	  if(sum(rags.ind)>0) # partial/completely missing
	  {
	    raggeds <- ragged[[seq(1,length(leads.rag))[rags.ind]]]
	    if(length(raggeds)<N) # partial missing
	    {
	      select.mat <- diag(N)[-raggeds,,drop=FALSE]
	      omit.mat <- diag(N)[raggeds,,drop=FALSE]
	    } else # completely missing
	    {
	      select.mat <- NULL
	      omit.mat <- diag(N)
	    }
	  }
	  non.raggeds <- setdiff(seq(1,N),raggeds)

	  #  get casts and covars of observations t.hash+1:t based on sigma-field_{t.hash+1:t}
	  # first, construct preds.x from known (unstored) entries and stored casts;
	  #   preds.x is E [ X_{t-1} | F_{t-1} ]
	  if(length(cast.index.t) > 0) { preds.x[,cast.index.t] <- casts.x }
	  # obtain E [ x_t | F_{t-1} ]
	  if(t.len > 0) { new.pred <- l.pred %*% matrix(preds.x[,(t-t.len):(t-1)],ncol=1) }

	  # second, get prediction variance, to obtain new.var given by Var [ x_t | F_{t-1} ]
	  #   cast.index.tlen tracks cast indices up to now, looking back t.len time points,
	  #   where t.len is the range at which all dependence is effectively nil
	  if(t.len > 0) { cast.index.tlen <- intersect(cast.indices,seq(t-t.len,t-1)) }
	  cast.len <- length(cast.index.tlen)
	  if(cast.len==0) # no casts within t.len time points
	  {
	    new.var <- v.pred
	  } else # at least one cast within t.len time points
	  {
	    casts.var.array <- array(casts.var,c(N,length(cast.index.t),N,length(cast.index.t)))
	    range.t <- (length(cast.index.t)-cast.len+1):length(cast.index.t)
	    casts.var.array <- casts.var.array[,,,range.t,drop=FALSE]
	    if(t-1-t.len>0) { l.pred <- cbind(matrix(0,N,N*(t-1-t.len)),l.pred) }
	    l.array <- array(l.pred,c(N,N,t-1))
	    l.array <- l.array[,,cast.index.tlen,drop=FALSE]
	    l.pred.tlen <- matrix(l.array,nrow=N)
	    new.var <- v.pred + l.pred.tlen %*%
	      matrix(casts.var.array[,range.t,,,drop=FALSE],nrow=cast.len*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	  }

	  # third, update casts.x by changing the stored portions of
	  #   E [ X_{t-1} | F_{t-1} ] to E [ X_{t-1} | F_t ] and
	  #   appending E [ x_t | F_t ] if partially/completely missing
	  if(cast.len>0)  # at least one cast within t.len time points
	  {
	    if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
	    {
	      new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	      update <- matrix(new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	                         (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred),nrow=N)
	      casts.x <- casts.x + update
	    }
	  }
	  if(length(raggeds)>0)   # add new cast E [ x_t | F_t ] if partially/completely missing
	  {
	    new.cast <- Re(z[,t,drop=FALSE])
	    partial.cast <- omit.mat %*% new.pred
	    if(length(raggeds)<N)  # partial missing case
	    {
	      partial.cast <- partial.cast + omit.mat %*% new.var %*% t(select.mat) %*%
	        solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
	    }
	    new.cast[raggeds,1] <- partial.cast
	    casts.x <- cbind(casts.x,new.cast)
	  }

	  # fourth, update casts.var by changing the stored portions of
	  #   Var [ X_{t-1} | F_{t-1} ] to Var [ X_{t-1} | F_t ] and
	  #   appending new covariances if partially/completely missing
	  if(cast.len>0)  # at least one cast within t.len time points
	  {
	    new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	    if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
	    {
	      update <- new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        select.mat %*% t(new.covar)
	      casts.var <- casts.var - update
	    }
	  }
	  if(length(raggeds)>0)   # add new covar [ x_t | F_t ] if partially/completely missing
	  {
	    proj <- diag(N) # completely missing case
	    if(length(raggeds)<N)  # partial missing case
	    {
	      proj <- diag(N) - t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        select.mat %*% new.var
	    }
	    # now augment casts.var
	    # special case: no cast within t.len time points.
	    #   Either: (i) this is the first cast,
	    #   Or: (ii) previous casts were long ago
	    if(cast.len==0) # special case
	    {
	      if(length(casts.var)==0)  # (i) of special case
	      {
	        new.block <- NULL
	      } else  # (ii) of special case
	      {
	        new.covar <- matrix(0,nrow=length(cast.index.t)*N,ncol=N)
	        new.block <- new.covar %*% proj
	      }
	    } else
	    {
	      new.block <- new.covar %*% proj
	    }
	    # do regular case (and special case) augmentation
	    casts.var <- rbind(cbind(casts.var,new.block),
	                       t(rbind(new.block,t(new.var %*% proj))))
	  }

	  # fifth, get ragged residuals
	  if(length(raggeds)>0)  # case of partially/completely missing
	  {
	    new.eps <- matrix(rep(1i,N),ncol=1)
	    new.det <- 1
	    if(length(raggeds)<N)  # partial missing case
	    {
	      partial.eps <- solve(t(chol(select.mat %*% new.var %*% t(select.mat)))) %*%
	        (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
	      new.eps[non.raggeds,1] <- partial.eps
	      new.det <- (det(chol(select.mat %*% new.var %*% t(select.mat))))^2
	    }
	  } else # case of full info
	  {
	    new.eps <- solve(t(chol(new.var))) %*% (Re(z[,t,drop=FALSE]) - new.pred)
	    new.det <- det(new.var)
	  }
	  eps <- rbind(eps,new.eps)
	  Qseq <- Qseq + t(Re(new.eps)) %*% Re(new.eps)
	  logdet <- logdet + log(new.det)

	  # updating
	  cast.index.t <- intersect(cast.indices,seq(t.hash+1,t))
	  preds.x <- Re(z[,1:t,drop=FALSE])

		#  get fore and aft predictors based on observations (t.hash+1):t
	  #  Notes: in stationary case, we obtain predictors for time t+1,
	  #   based on data over times t.hash+1,...t.
	  #   That is, the predictors computed below are used in the next
	  #   increment of the t-loop.
	  if(t==(t.hash+1))
		{

			u.seq <- solve(x.acf[,1,]) %*% x.acf[,2,]
			l.seq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			gam.Seq <- x.acf[,2,]
			gam.Flip <- x.acf[,2,]
			c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
			d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
			a.seq <- u.seq
			b.seq <- l.seq

			u.qes <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			l.qes <- solve(x.acf[,1,]) %*% x.acf[,2,]
			gam.qeS <- t(x.acf[,2,])
			gam.pilF <- t(x.acf[,2,])
			c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
			d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
			a.qes <- u.qes
			b.qes <- l.qes

		} else # case of t.hash+2 <= t <= T
		{
			if((t.star == T) && (t<T))
			{

				pacf <- x.acf[,t+1-t.hash,] - gam.Seq %*% u.seq
				l.factor <- solve(c.mat) %*% t(pacf)
				new.l <- l.seq - u.seq %*% l.factor
				u.factor <- solve(d.mat) %*% pacf
				new.u <- u.seq - l.seq %*% u.factor
				l.seq <- rbind(l.factor,new.l)
				u.seq <- rbind(new.u,u.factor)
				gam.Seq <- cbind(x.acf[,t+1-t.hash,],gam.Seq)
				gam.Flip <- rbind(gam.Flip,x.acf[,t+1-t.hash,])
				c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
				d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
				a.next <- a.seq - b.seq %*% u.factor
				b.next <- b.seq - a.seq %*% l.factor
				a.seq <- rbind(a.next,0*diag(N)) + (c(rep(0,t-1-t.hash),1) %x% diag(N)) %*% u.factor
				b.seq <- rbind(0*diag(N),b.next) + (c(1,rep(0,t-1-t.hash)) %x% diag(N)) %*% l.factor

				fcap <- t(x.acf[,t+1-t.hash,]) - gam.qeS %*% u.qes
				l.rotcaf <- solve(c.tam) %*% t(fcap)
				wen.l <- l.qes - u.qes %*% l.rotcaf
				u.rotcaf <- solve(d.tam) %*% fcap
				wen.u <- u.qes - l.qes %*% u.rotcaf
				l.qes <- rbind(l.rotcaf,wen.l)
				u.qes <- rbind(wen.u,u.rotcaf)
				gam.qeS <- cbind(t(x.acf[,t+1-t.hash,]),gam.qeS)
				gam.pilF <- rbind(gam.pilF,t(x.acf[,t+1-t.hash,]))
				c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
				d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
				a.next <- a.qes - b.qes %*% u.rotcaf
				b.next <- b.qes - a.qes %*% l.rotcaf
				a.qes <- rbind(0*diag(N),a.next) + (c(1,rep(0,t-1-t.hash)) %x% diag(N)) %*% u.rotcaf
				b.qes <- rbind(b.next,0*diag(N)) + (c(rep(0,t-1-t.hash),1) %x% diag(N)) %*% l.rotcaf

			}
			if(sqrt(sum(diag(l.factor %*% t(l.factor)))) < thresh)
			{
				t.star <- min(t.star,t)
			}
		}
		t.len <- dim(b.seq)[1]/N
		l.pred <- t(b.seq)
		l.derp <- t(b.qes)
		v.pred <- (x.acf[,1,] - matrix(aperm(x.acf[,(t.len+1):2,,drop=FALSE],c(1,3,2)),nrow=N) %*% l.seq)
		v.derp <- (x.acf[,1,] - matrix(aperm(x.acf[,(t.len+1):2,,drop=FALSE],c(3,1,2)),nrow=N) %*% l.qes)
	} }

	# Backward Pass:
	if(t.hash > 0) {
	for(t in (t.hash):1)
	{
	  # determine whether full info, or partial/completely missing
	  #  base case: full info
	  select.mat <- diag(N)
	  omit.mat <- NULL
	  raggeds <- NULL
	  rags.ind <- leads.rag %in% t
	  if(sum(rags.ind)>0) # partial/completely missing
	  {
	    raggeds <- ragged[[seq(1,length(leads.rag))[rags.ind]]]
	    if(length(raggeds)<N) # partial missing
	    {
	      select.mat <- diag(N)[-raggeds,,drop=FALSE]
	      omit.mat <- diag(N)[raggeds,,drop=FALSE]
	    } else # completely missing
	    {
	      select.mat <- NULL
	      omit.mat <- diag(N)
	    }
	  }
	  non.raggeds <- setdiff(seq(1,N),raggeds)

	  #  get casts and covars of observations t:T based on sigma-field_{t:T}
	  # first, construct preds.x from known (unstored) entries and stored casts;
	  #   preds.x is E [ X_{t-1} | F_{t+1} ]
	  if(length(cast.index.t) > 0) { preds.x[,cast.index.t] <- casts.x }
	  # obtain E [ x_t | F_{t-1} ]
	  if(t.len > 0) { new.pred <- l.derp %*% matrix(preds.x[,(t+1):(t+t.len)],ncol=1) }

	  # second, get prediction variance, to obtain new.var given by Var [ x_t | F_{t+1} ]
	  #   cast.index.tlen tracks cast indices up to now, looking forward t.len time points,
	  #   where t.len is the range at which all dependence is effectively nil
	  if(t.len > 0) { cast.index.tlen <- intersect(cast.indices,seq(t+1,t+t.len)) }
	  cast.len <- length(cast.index.tlen)
	  if(cast.len==0) # no casts within t.len time points
	  {
	    new.var <- v.derp
	  } else # at least one cast within t.len time points
	  {
	    casts.var.array <- array(casts.var,c(N,length(cast.index.t),N,length(cast.index.t)))
#	    range.t <- (length(cast.index.t)-cast.len+1):length(cast.index.t)
	    range.t <- 1:cast.len
	    casts.var.array <- casts.var.array[,,,range.t,drop=FALSE]
	    l.array <- array(l.derp,c(N,N,T-t))
	    l.array <- l.array[,,cast.index.tlen-t,drop=FALSE]
	    l.pred.tlen <- matrix(l.array,nrow=N)
	    new.var <- v.derp + l.pred.tlen %*%
	      matrix(casts.var.array[,range.t,,,drop=FALSE],nrow=cast.len*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	  }

	  # third, update casts.x by changing the stored portions of
	  #   E [ X_{t+1} | F_{t+1} ] to E [ X_{t+1} | F_t ] and
	  #   appending E [ x_t | F_t ] if partially/completely missing
	  if(cast.len>0)  # at least one cast within t.len time points
	  {
	    if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
	    {
	      new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	      update <- matrix(new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	                         (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred),nrow=N)
	      casts.x <- casts.x + update
	    }
	  }
	  if(length(raggeds)>0)   # add new cast E [ x_t | F_t ] if partially/completely missing
	  {
	    new.cast <- Re(z[,t,drop=FALSE])
	    partial.cast <- omit.mat %*% new.pred
	    if(length(raggeds)<N)  # partial missing case
	    {
	      partial.cast <- partial.cast + omit.mat %*% new.var %*% t(select.mat) %*%
	        solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
	    }
	    new.cast[raggeds,1] <- partial.cast
	    casts.x <- cbind(new.cast,casts.x)
	  }

	  # fourth, update casts.var by changing the stored portions of
	  #   Var [ X_{t+1} | F_{t+1} ] to Var [ X_{t+1} | F_t ] and
	  #   appending new covariances if partially/completely missing
	  if(cast.len>0)  # at least one cast within t.len time points
	  {
	    new.covar <- matrix(casts.var.array,nrow=length(cast.index.t)*N,ncol=cast.len*N) %*% t(l.pred.tlen)
	    if(length(raggeds)<N)  # update only if full info or partially missing (do nothing if fully missing)
	    {
	      update <- new.covar %*% t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        select.mat %*% t(new.covar)
	      casts.var <- casts.var - update
	    }
	  }
	  if(length(raggeds)>0)   # add new covar [ x_t | F_t ] if partially/completely missing
	  {
	    proj <- diag(N) # completely missing case
	    if(length(raggeds)<N)  # partial missing case
	    {
	      proj <- diag(N) - t(select.mat) %*% solve(select.mat %*% new.var %*% t(select.mat)) %*%
	        select.mat %*% new.var
	    }
	    # now augment casts.var
	    # special case: no cast within t.len time points.
	    #   Either: (i) this is the first cast,
	    #   Or: (ii) previous casts were long ago
	    if(cast.len==0) # special case
	    {
	      if(length(casts.var)==0)  # (i) of special case
	      {
	        new.block <- NULL
	      } else  # (ii) of special case
	      {
	        new.covar <- matrix(0,nrow=length(cast.index.t)*N,ncol=N)
	        new.block <- new.covar %*% proj
	      }
	    } else
	    {
	      new.block <- new.covar %*% proj
	    }
	    # do regular case (and special case) augmentation
	    casts.var <- rbind(t(rbind(t(new.var %*% proj),new.block)),
	                       cbind(new.block,casts.var))
	  }

	  # fifth, get ragged residuals
	  if(length(raggeds)>0)  # case of partially/completely missing
	  {
	    new.eps <- matrix(rep(1i,N),ncol=1)
	    new.det <- 1
	    if(length(raggeds)<N)  # partial missing case
	    {
	      partial.eps <- solve(t(chol(select.mat %*% new.var %*% t(select.mat)))) %*%
	        (Re(z[non.raggeds,t,drop=FALSE]) - select.mat %*% new.pred)
	      new.eps[non.raggeds,1] <- partial.eps
	      new.det <- (det(chol(select.mat %*% new.var %*% t(select.mat))))^2
	    }
	  } else # case of full info
	  {
	    new.eps <- solve(t(chol(new.var))) %*% (Re(z[,t,drop=FALSE]) - new.pred)
	    new.det <- det(new.var)
	  }
	  eps <- rbind(eps,new.eps)
	  Qseq <- Qseq + t(Re(new.eps)) %*% Re(new.eps)
	  logdet <- logdet + log(new.det)

	  # updating
	  cast.index.t <- intersect(cast.indices,seq(t,T))
	  preds.x <- Re(z[,1:T,drop=FALSE])

		#  get predictors based on observations t+1:T
	  #  Notes: in stationary case, we obtain predictors for time t-1,
	  #   based on data over times t,...T.
	  #   That is, the predictors computed below are used in the next
	  #   increment of the t-loop.
	  if(t==T)	# if no forward pass happened, initialize
		{

			u.seq <- solve(x.acf[,1,]) %*% x.acf[,2,]
			l.seq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			gam.Seq <- x.acf[,2,]
			gam.Flip <- x.acf[,2,]
			c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
			d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
			a.seq <- u.seq
			b.seq <- l.seq

			u.qes <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
			l.qes <- solve(x.acf[,1,]) %*% x.acf[,2,]
			gam.qeS <- t(x.acf[,2,])
			gam.pilF <- t(x.acf[,2,])
			c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
			d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
			a.qes <- u.qes
			b.qes <- l.qes

		} else # case of 1 <= t <= t.hash and t.hash < T
		{
			if((t.star == T) && (t>1))
			{

				pacf <- x.acf[,T-t+2,] - gam.Seq %*% u.seq
				l.factor <- solve(c.mat) %*% t(pacf)
				new.l <- l.seq - u.seq %*% l.factor
				u.factor <- solve(d.mat) %*% pacf
				new.u <- u.seq - l.seq %*% u.factor
				l.seq <- rbind(l.factor,new.l)
				u.seq <- rbind(new.u,u.factor)
				gam.Seq <- cbind(x.acf[,T-t+2,],gam.Seq)
				gam.Flip <- rbind(gam.Flip,x.acf[,T-t+2,])
				c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
				d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
				a.next <- a.seq - b.seq %*% u.factor
				b.next <- b.seq - a.seq %*% l.factor
				a.seq <- rbind(a.next,0*diag(N)) + (c(rep(0,T-t),1) %x% diag(N)) %*% u.factor
				b.seq <- rbind(0*diag(N),b.next) + (c(1,rep(0,T-t)) %x% diag(N)) %*% l.factor

				fcap <- t(x.acf[,T-t+2,]) - gam.qeS %*% u.qes
				l.rotcaf <- solve(c.tam) %*% t(fcap)
				wen.l <- l.qes - u.qes %*% l.rotcaf
				u.rotcaf <- solve(d.tam) %*% fcap
				wen.u <- u.qes - l.qes %*% u.rotcaf
				l.qes <- rbind(l.rotcaf,wen.l)
				u.qes <- rbind(wen.u,u.rotcaf)
				gam.qeS <- cbind(t(x.acf[,T-t+2,]),gam.qeS)
				gam.pilF <- rbind(gam.pilF,t(x.acf[,T-t+2,]))
				c.tam <- x.acf[,1,] - t(gam.pilF) %*% u.qes
				d.tam <- x.acf[,1,] - gam.qeS %*% l.qes
				a.next <- a.qes - b.qes %*% u.rotcaf
				b.next <- b.qes - a.qes %*% l.rotcaf
				a.qes <- rbind(0*diag(N),a.next) + (c(1,rep(0,T-t)) %x% diag(N)) %*% u.rotcaf
				b.qes <- rbind(b.next,0*diag(N)) + (c(rep(0,T-t),1) %x% diag(N)) %*% l.rotcaf

			}
			if(sqrt(sum(diag(l.factor %*% t(l.factor)))) < thresh)
			{
				t.star <- min(t.star,t)
			}
		}
		t.len <- dim(b.seq)[1]/N
		l.pred <- t(b.seq)
		l.derp <- t(b.qes)
		v.pred <- (x.acf[,1,] - matrix(aperm(x.acf[,(t.len+1):2,,drop=FALSE],c(1,3,2)),nrow=N) %*% l.seq)
		v.derp <- (x.acf[,1,] - matrix(aperm(x.acf[,(t.len+1):2,,drop=FALSE],c(3,1,2)),nrow=N) %*% l.qes)
	} }

	}

	eps <- matrix(eps,nrow=N)
	if(length(cast.index.t) > 0) { preds.x[,cast.index.t] <- casts.x }
	lik <- Qseq + logdet

	if(debug) { print(lik) }

	return(list(casts.x,casts.var,c(Qseq,logdet),eps))
}

