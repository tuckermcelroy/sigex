mvar.forecast <- function(x.acf,z,needMSE)
{

	##########################################################################
	#
	#	mvar.forecast
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
	#	Purpose: compute multi-step forecasts and predictors of a multivariate process
	#		via Levinson-Durbin algorithm
	#	Background:	
	#		A multivariate difference-stationary process x_t with
	#			w_t = delta(B) x_t 
	#		may be observed with missing values, and one wants to compute
	#		Gaussian conditional expectations of missing values (midcasts),
	#		or future values (forecasts), or past values (aftcasts).
	#		Also of interest is the Gaussian likelihood resulting from
	#		such a sample, and the residuals.  
	#	Inputs:
	#		x.acf: array of dimension N x T x N of autocovariances for process w_t,
	#			where there are N series, of total length T each.
	#		z: differenced data as N x (T+H) matrix, with missing values at
	#			various time points.  Presumes first T observations are not missing, 
	#			and latter H observations are missing, being encoded
	#			with 1i in that entry.  That is, 
	#			Im(z[,t]) = rep(1i,N) encodes missing values.
	#		needMSE: a binary flag, with 1 indicating that forecast MSE should be computed;
	#			else with a 0 this will be skipped (runs faster in this mode).
 	#	Outputs:
	#		list containing preds.x and pred.stack
	#		preds.x: N x (T+H) matrix of data with forecasts, where H
	#			is the total number of time indices with missing values.
	#		pred.stack: NT x H matrix of predictors, which can be used to
	#			obtain uncertainty...
	#	Notes: running this code with needMSE=1 makes it slower, but yields
	#		pred.stack, however the routine mvar.midcast is preferred for
	#		obtaining any casts with uncertainty.  I think the pred.stack feature
	#		here should be deprecated, but not sure if I want to completely
 	#		remove it in case another application comes up.
	#
	####################################################################
		
	thresh <- 10^(-16)
	N <- dim(z)[1]
	TH <- dim(z)[2]
	all.series <- seq(1,N)
	all.indices <- seq(1,TH)
	full.indices <- all.indices[colSums(z==1i)==0]
	cast.indices <- setdiff(all.indices,full.indices)
	H <- length(cast.indices)
	T <- length(full.indices)

	# initialization
	preds.x <- Re(z[,1:T,drop=FALSE])
	pred.stack <- NULL

	# first pass updates
	u.seq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	l.seq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gam.Seq <- x.acf[,2,]
	gam.Flip <- x.acf[,2,]
	c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
	d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
	
	# looping
	t.star <- T
	for(t in 1:(TH-2))
	{
		if(t.star == T)
		{
			pacf <- x.acf[,t+2,] - gam.Seq %*% u.seq
			l.factor <- solve(c.mat) %*% t(pacf)
			new.l <- l.seq - u.seq %*% l.factor
			u.factor <- solve(d.mat) %*% pacf
			new.u <- u.seq - l.seq %*% u.factor
			l.seq <- rbind(l.factor,new.l)
			u.seq <- rbind(new.u,u.factor)
			gam.Seq <- cbind(x.acf[,t+2,],gam.Seq)
			gam.Flip <- rbind(gam.Flip,x.acf[,t+2,])
			c.mat <- x.acf[,1,] - t(gam.Flip) %*% u.seq
			d.mat <- x.acf[,1,] - gam.Seq %*% l.seq
		}
		if((sqrt(sum(diag(l.factor %*% t(l.factor)))) < thresh) && (t < T-1)) 
		{ 
			t.star <- min(t.star,t+1)
		}
		if(t==(T-1))
		{
			if(needMSE) {
				pred.next <- l.seq
				pred.stack <- pred.next
			} else {
				pred.next <- l.seq
				t.len <- dim(l.seq)[1]/N
				pred.x <- t(l.seq) %*% matrix(Re(z[,(t+1-t.len+1):(t+1)]),ncol=1)
				preds.x <- cbind(preds.x,pred.x)
			}
		}
		if(t > (T-1)) 
		{
			if(needMSE) {
				t.len <- dim(l.seq)[1]/N
				stack.dim <- dim(pred.stack)
				pred.next <- cbind(diag(stack.dim[1]),pred.stack)[,(stack.dim[1]-t.len+stack.dim[2]+1):(stack.dim[1]+stack.dim[2]),drop=FALSE] %*% l.seq
				pred.stack <- cbind(pred.stack,pred.next)
			} else {
				pred.next <- l.seq
				t.len <- dim(l.seq)[1]/N
				pred.x <- t(l.seq) %*% matrix(Re(preds.x[,(t+1-t.len+1):(t+1)]),ncol=1)
				preds.x <- cbind(preds.x,pred.x)
			}
		}
	}

	if(needMSE) {
		x.cast <- t(pred.stack) %*% matrix(Re(z[,(T-t.star+1):T]),ncol=1)
		x.cast <- matrix(x.cast,nrow=N)
		preds.x <- cbind(Re(z[,1:T,drop=FALSE]),x.cast)
	} 
	
	return(list(preds.x,pred.stack)) 
}
