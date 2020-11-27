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
IntegerVector getRagIndex(ComplexVector zval)
{
  int N = zval.length();
  arma::vec imPart = getIm(zval);
  IntegerVector rag_series;
  for(int k = 0; k < N; k++)
  {
    if(imPart(k) == 1)  // This means zval[k] is a missing value
    { 
      rag_series.push_back(k);  
    }
  }
  return rag_series;
}

// [[Rcpp::export]]
IntegerVector subsetting(IntegerVector x, IntegerVector y, double z) 
{
  return x[y > z];
}

// [[Rcpp::export]] 
arma::mat matmult(arma::mat A, arma::mat B) { return A % B; }

// [[Rcpp::export]]
arma::mat kronprod(arma::mat A, arma::mat B) { return kron(A,B); }

// [[Rcpp::export]] 
arma::vec vec2mat(arma::vec x) { return(x) ; }

// [[Rcpp::export]]
List mvar_midcast(arma::cube x_acf, ComplexMatrix z, NumericVector delta, bool debug) 
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
//		x_acf: array of dimension N x N x T of autocovariances for process w_t,
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

//  double thresh = 10E-16;
  int N = z.nrow();
  int T = z.ncol();
  arma::mat idmat = arma::eye<arma::mat>(N,N);

  IntegerVector all_series (N,1);
  for(int i = 0; i < N; i++) { all_series [i] = i+1; } 
  // All indices for the time series' extended sample
  IntegerVector all_indices (T,1);
  // Subset of indices for which there is no missing values
  IntegerVector full_indices;
  // Subset of indices for which there is at least one missing value
  IntegerVector cast_indices;
  // ragged is N x T  boolean matrix, with FALSE if corresponding
  //  value in z is missing.  In the R version, ragged was a list
  //  object but cpp can't do it.  
  LogicalMatrix ragged (N,T);
  ragged.fill(TRUE);
  
  // This loop computes full_indices and cast_indices from all_indices,
  //  by determining whether there are any missing values at index t.
  //  Note: cast_indices and full_indices can be empty (length zero vector).
  //  rag_series gets component indices of missing values at index t,
  //   by checking whether the imaginary part of column t of z equals 1.
  for(int t = 0; t < T; t++) 
  { 
    int tp1 = t+1;
    all_indices[t] = tp1; 
    IntegerVector rag_series = getRagIndex(z(_,t));
    if(rag_series.length() == 0)  // In this case, z[_,t] had no missing values
    { 
      full_indices.push_back(tp1);
    } else              // In this case, z[,t] had some missing values
    {
      cast_indices.push_back(tp1);
      for(int k = 0; k < rag_series.length(); k++)
      {
        int j = rag_series[k];
        ragged(j,t) = FALSE;
      }
    }
  }
  
  // d is the degree of the differencing polynomial
  int d = delta.length() - 1;
 
// This version does not presume that the first d values are not missing 
// Find t.hash, earliest time for which d contiguous values follow,
//	such that data at times t.hash+1,...,t.hash+d exists
  int t_hash = 0;
  IntegerVector ind_data = cast_indices;
  ind_data.push_front(0);
  ind_data.push_back(T+1);
  IntegerVector gaps = diff(ind_data);
  gaps.push_back(0);
  if(max(gaps) > d) 
  { 
     t_hash = min(subsetting(ind_data,gaps,d));
  } else
  {
    Rcerr << "There are not " << d << " contiguous values\n";
  }  

//  Rcout << "t.hash  "  << t_hash;
  
//  Main code block for nonstationary case (d > 0)      
  if(d > 0) 
  {

// t = t_hash + d case as initialization
//  get predictors based on observations t_hash+1 to t
    arma::mat l_pred (N,N*d);
    arma::vec del_aft = -1*rev(tail(delta,d))/delta(0);
    arma::mat del_aft_mat = vec2mat(del_aft);
    l_pred = kronprod(del_aft_mat,idmat);
    l_pred = l_pred.t();
    arma::mat l_derp (N,N*d);
    arma::vec del_fore = -1*rev(head(delta,d))/delta(d);
    arma::mat del_fore_mat = vec2mat(del_fore);
    l_derp = kronprod(del_fore_mat,idmat);
    l_derp = l_derp.t();
    arma::mat v_pred = x_acf.slice(0)/(delta(0)*delta(0));
    arma::mat v_derp = x_acf.slice(0)/(delta(d)*delta(d));
    
//  get casts and covars of observations t_hash+1 to t based
//    on sigma-field_{t_hash+1:t}
//   Note: store more than needed in preds_x, makes it easier 
//    for indexing later
    arma::mat preds_x (N,t_hash+d);
    for(int t = 0; t < t_hash+d; t++)
    {
      preds_x.col(t) = getRe(z(_,t));
    }
    arma::mat new_covar;
    arma::mat casts_x;
    arma::mat casts_var;
    arma::mat eps;
    double Qseq = 0;
    double logdet = 0;

// track indices of casted variables up to present time
    IntegerVector init_index = seq(t_hash+1,t_hash+d);
    IntegerVector cast_index_t = intersect(cast_indices,init_index);
    int t_star = T;
    int t_len = d;
    
// Forward Pass:	
    if(t_hash < T-d) 
    {
      for(int t = t_hash+d+1; t <= T; t++)
      { 
        int tm1 = t-1;
// determine whether full info, or partial/completely missing
//  base case: full info
        arma::mat select_mat;
        arma::mat omit_mat;
        LogicalVector rag_vec = ragged(_,tm1);          
        IntegerVector raggeds;
        IntegerVector non_raggeds;
        for(int k = 0; k < N; k++)
        {
          if(rag_vec(k)) { non_raggeds.push_back(k+1); } else { raggeds.push_back(k+1); }
        }  
//        if(raggeds.length() > 0) { omit_mat = idmat.submat(raggeds,seq(1,N)); }
//        if(non_raggeds.length() > 0) { select_mat = idmat.submat(non_raggeds,seq(1,N)); }
//if(raggeds.length() > 0) { omit_mat = idmat.submat(1,2,1,1); }

          
    Rcout << "select matrix "  << select_mat << "\n";
      Rcout << "omit matrix "  << omit_mat << "\n";
      Rcout << "raggeds " << raggeds << "\n";       
    Rcout << "non raggeds " << non_raggeds << "\n";       
    
        
      } 
    }  // end forward pass clause
      
  }   // end d > 0  clause
   
  List L = List::create(ragged);
  
  return L;

}