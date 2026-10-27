// Use the RcppArmadillo package 
// Requires different header file from Rcpp.h 
#include <RcppArmadillo.h> 
#include <algorithm> // Required for std::set_difference
#include <iterator>  // Required for std::back_inserter
#include <complex>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getRe(ComplexVector x)
{
  arma::cx_vec x_arma = as<arma::cx_vec>(x);
  arma::vec x_real = arma::real(x_arma);
  return x_real;
}

// [[Rcpp::export]]
arma::vec getIm(ComplexVector x) 
{
  arma::cx_vec x_arma = as<arma::cx_vec>(x);
  arma::vec x_imag = arma::imag(x_arma);
  return x_imag;  
}

// [[Rcpp::export]]
// Returns the components using C convention with imaginary part == 1
//   This is used as a flag for missing value 
IntegerVector get_rag_index(ComplexVector zval)
{
  int N = zval.length();
  arma::vec imPart = getIm(zval);
  IntegerVector rag_series;
  for(int k = 0; k < N; k++)
  {
    if(imPart(k) == 1)  // This means element k of zval is a missing value
    { 
      rag_series.push_back(k);  
    }
  } // else no missing values in zval, and rag_series is empty
  return rag_series;
}

// [[Rcpp::export]] 
// Computes Hadamard product of matrices A and B of the same dimension
arma::mat matmult(arma::mat A, arma::mat B) { return A % B; }

// [[Rcpp::export]]
// Computes kronecker i.e. tensor product of matrices A and B
arma::mat kronprod(arma::mat A, arma::mat B) { return kron(A,B); }

// [[Rcpp::export]] 
// Formats a vector as a matrix with column dimension == 1
arma::vec vec2mat(arma::vec x) { return(x) ; }

// [[Rcpp::export]]
int mvar_midcast(arma::cube x_acf, ComplexMatrix z, NumericVector delta, bool debug = false) 
{
  
  //***********************************************************************
  //
  //	mvar_midcast.cpp
  // 	    Copyright (C) 2026  Tucker McElroy
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
  int thresh = -1;
  int N = z.nrow();
  int T = z.ncol();
  arma::mat idmat = arma::eye<arma::mat>(N,N);
  
  IntegerVector all_series (N,1);
  for(int i = 0; i < N; i++) { all_series [i] = i+1; } 
  
  // All indices for the time series' extended sample
  IntegerVector all_indices (T,1);
  for(int i = 0; i < T; i++) { all_indices [i] = i+1; } 
  
  // Subset of indices for which there are no missing values
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
    IntegerVector rag_series = get_rag_index(z(_,t));
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
  NumericVector delta_rev = rev(delta);
  arma::vec delta_flip = as<arma::vec>(delta_rev);
  arma::mat delta_lead = kronprod(delta_flip,idmat);
  
  // Get real part of the z matrix
  arma::cx_mat z_arma = as<arma::cx_mat>(z);
  arma::mat z_real = arma::real(z_arma);

  // We do not presume that the first d values are not missing 
  // Find t.hash, earliest time for which d contiguous values follow,
  //	such that data at times t.hash+1,...,t.hash+d exists
  int t_hash = 0;
  IntegerVector ind_data = cast_indices;
  ind_data.push_front(0);
  ind_data.push_back(T+1);
  IntegerVector gaps = diff(ind_data);
  gaps.push_back(0);
  
  int gaps_max = max(gaps);
  if(gaps_max > d)   
  {
    // write loop, need C indexing to work with R indexing!
    NumericVector filtered = ind_data[gaps > d];
    
//    t.hash <- min(ind.data[gaps > d]) 
  }  
    
  return gaps_max;
  
}  

