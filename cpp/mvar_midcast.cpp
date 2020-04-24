// Use the RcppArmadillo package 
// Requires different header file from Rcpp.h 
#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getRe(arma::cx_vec x) {
  return arma::real(x);
}

// [[Rcpp::export]]
arma::vec getIm(arma::cx_vec x) {
  return arma::imag(x);
}

// [[Rcpp::export]]
List mvar_midcast(List x_acf, ComplexMatrix z, NumericVector delta, bool debug) 
{
  
//***********************************************************************
//
//	mvar.midcast.cpp
// 	    Copyright (C) 2020  Tucker McElroy
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//**************************************************************************

//************************ Documentation **********************************
//
//	Purpose: compute multi-step imputations and predictors of a 
//    multivariate process
//		via Levinson-Durbin algorithm with  missing values
//	Background:	
//		A multivariate difference-stationary process x_t with
//			w_t = delta(B) x_t 
//		may be observed with missing values, and one wants to compute
//		Gaussian conditional expectations of missing values (midcasts),
//		or future values (forecasts), or past values (aftcasts).
//		Also of interest is the Gaussian likelihood resulting from
//		such a sample, and the residuals.  
//		It is required that at least d
//		contiguous values be observed, where d is the order of delta(B).
//		If the first d values are contiguous, we can do a forward pass;
//		otherwise, the first set of d contiguous values starts after
//		time index 1, and is marked by t.hash.  In this case, we must
//		also do a backward pass, involving aftcasts.
//	Inputs:
//		x.acf: array of dimension N x T x N of autocovariances for process w_t,
//			where there are N series, of total length T each.
//		z: raw data as N x T matrix with missing values at various time points.
//			Missing values are at any 1 <= t <= T, and occur for some of the N series
//     (the ragged case), and are denoted with a 1i.  That is, 
//			Im(z[,t]) = rep(1i,N) or subset thereof encodes missing values.
//		delta: differencing polynomial (corresponds to delta(B) in Background)
//			written in format c(delta0,delta1,...,deltad)
//   debug: set to TRUE if lik values should be printed to screen
//	Notes: to get H forecasts, append matrix(1i,N,H) to input x.  To get aftcasts,
//		prepend the same.  T will be the second dimension of z, and includes
//		the spots taken by aftcasts and forecasts.  (So the T for the original
//		dataset could be less than the T used in this function.)
//	Outputs:
//		list containing casts.x, casts.var, c(Qseq,logdet), and eps
//		casts.x: N x H matrix of backcasts, midcasts, aftcasts, where H
//			is the total number of time indices with missing values.
//			So if times.na is a subset of seq(1,T) corresponding to indices
//			with missing values, we can fill in all NAs via	z[,times.na] <- casts.x.
//     If a subset is missing at some time t, casts.x for that time t contains
//     the casts together with the known values.
//		casts.var: NH x NH matrix of covariances of casting errors.
//			note that casts.var.array <- array(casts.var,c(N,H,N,H)) 
//			corresponds to cast.var.array[,j,,k] equal to the 
//			covariance between the jth and kth casting errors.
//     If a subset is missing at some time t, the corresponding block of casts.var
//     will have zeros corresponding to the known values.
//		Qseq: quadratic form portion of the Gaussian divergence based on
//		missing value formulation (McElroy and Monsell, 2015 JASA)
//		logdet: log determinant portion of the Gaussian divergence
//		eps: residuals from casting recursions, defined in the manner
//			discussed in Casting paper.  Dimension is N x T-(H+d)
//
//####################################################################

  double thresh = 10E-16;
  int N = z.nrow();
  int T = z.ncol();
   
  IntegerVector all_series (N,1);
  for(int i = 0; i < N; i++) { all_series [i] = i+1; } 
  IntegerVector all_indices (T,1);
  LogicalVector full_count (T,1);
  LogicalVector cast_count (T,1);
  ComplexVector zval;
  arma::vec imPart;
  bool is_na;
  
  for(int t = 0; t < T; t++) 
  { 
    all_indices[t] = t+1; 
    zval = z(_,t);
    imPart = getIm(zval);
    is_na = FALSE;
    full_count(t) = FALSE;
    cast_count(t) = TRUE;
    for(int k = 0; k < N; k++)
    {
      if(imPart(k) == 1) { is_na = TRUE; }
    }
    if(is_na == FALSE) 
    { 
      full_count(t) = TRUE; 
      cast_count(t) = FALSE;
    }
  }
  
  Rcout << full_count << "\n";
  IntegerVector full_indices;
  
  // if(any(full_count) == TRUE) 
  // { 
  //   full_indices = ifelse(full_count,all_indices,0); 
  // } 
  // 
  // Rcout << full_indices << "\n";
  
//  ragged <- list()
//  leads.rag <- NULL
//  for(t in 1:length(cast.indices))
//  {
//   rag.series <- all.series[Im(z[,cast.indices[t]])==1,drop=FALSE]
//    if(length(rag.series)<=N) 
//    { 
//      ragged[[length(ragged)+1]] <- rag.series 
//      leads.rag <- c(leads.rag,cast.indices[t])
//    }
//  }
//  d <- length(delta) - 1
  
  List L = List::create(all_series, all_indices, thresh);
  
  return L;

}