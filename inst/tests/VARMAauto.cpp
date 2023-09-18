// Use the RcppArmadillo package 
// Requires different header file from Rcpp.h 
#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube polymul_mat(arma::cube amat,arma::cube bmat)
{

  int p = amat.n_slices -1;
  int q = bmat.n_slices -1;
  int N = amat.n_cols;
  
  int r = p+q;
  arma::cube bmat_pad(N,N,r+1);
  arma::cube cmat(N,N,r+1);
  
  for(int i = 0; i < q+1; i++) 
  { 
    bmat_pad.slice(i) = bmat.slice(i); 
  }
  cmat.slice(0) = amat.slice(0) * bmat_pad.slice(0);
  
  if(r > 0) {
  for(int j = 1; j < r+1; j++)
  {
    cmat.slice(j) = amat.slice(0) * bmat_pad.slice(j);
    if(p > 0) 
    {
      int l = std::min(p,j);
      for(int k = 0; k < l; k++)
      { 
        arma::mat temp = amat.slice(k+1) * bmat_pad.slice(j-k-1);
        cmat.slice(j) = cmat.slice(j) + temp;
      }
    }  
  }}
  
  return cmat;
  
}
 

// [[Rcpp::export]]
arma::cube VARMA_auto(arma::mat param,int p, int q, int maxlag)
{
  
  //***********************************************************************
  //
  //	VARMAauto.cpp
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
  //	Purpose: computes autocovariances of VARMA
  //	Background: function computes autocovariances of VARMA (p,q) from lag zero
  //		to maxlag, with array inputs phi and theta.  VARMA equation:
  //  	(1 - phi[1]B ... - phi[p]B^p) X_t = (1 + theta[1]B ...+ theta[q]B^q) WN_t
  //	Note: for absent VAR or VMA portions, pass in NULL
  //	Inputs:
  //    param: matrix of dimension m x (p+q+1)m, equals [ phi | theta | sigma ]
  //		  phi: block matrix of dimension m x mp of VAR coefficients
  //		  theta: block matrix of dimension m x mq of VMA coefficients 
  //		  sigma: m x m covariance matrix of white noise
  //  Notes: for pure VMA, leave out phi; for pure VAR, leave out theta.
  //	Outputs:
  //		autocovariances at lags 0 through maxlag, as array of dimension m x m x (maxlag+1)
  //
  //**************************************************************************
  
  int m = param.n_rows;
  arma::mat idmat = arma::eye<arma::mat>(m,m);
  arma::mat phi = idmat;  // only to get phi globally defined, but if p = 0, phi should be NULL
  arma::cube gam_ma(m,m,q+1);
  arma::mat sigma = param.cols((p+q)*m,(p+q+1)*m-1);
  arma::cube sigma_cube(m,m,1);
  sigma_cube.slice(0) = sigma;
  arma::cube gamma_final(m,m,maxlag+1);
  
  arma::mat Kmat (m*m,m*m);
  Kmat.fill(0);
  for(int j = 0; j < m; j++)
  {
    for(int k = 0; k < m; k++)
    {
      int row_index = m*j + k;
      int col_index = m*k + j;
      Kmat(row_index,col_index) = 1;
    }
  }
  
  if(q == 0) 
  { 
    gam_ma.slice(0) = sigma;
  } else
  {
    arma::mat theta = param.cols(p*m,(p+q)*m-1);
    arma::cube theta_cube(m,m,q+1);
    theta_cube.slice(0) = idmat;
    for(int i = 0; i < q; i++)
    {
      theta_cube.slice(i+1) = theta.cols(i*m,(i+1)*m-1);
    }
    arma::cube prod_cube = polymul_mat(theta_cube,sigma_cube);
    theta_cube.slice(q) = idmat;
    for(int i = 0; i < q; i++)
    {
      arma::mat temp = theta.cols(i*m,(i+1)*m-1);
      theta_cube.slice(q-1-i) = temp.t();
    }
    arma::cube out = polymul_mat(prod_cube,theta_cube);
    gam_ma = out.slices(q,2*q);
  }
  arma::vec gamvec_ma = vectorise(gam_ma);
  arma::cube gam_mixcube(m,m,q+1);
  arma::cube gam_armacube(m,m,p+1);
  
  if(p > 0)
  {
    arma::mat id2mat = arma::eye<arma::mat>(m*m,m*m);
    phi = param.cols(0,p*m-1);
    arma::mat Amat(m*m*(p+1),m*m*(2*p+1));
    Amat.fill(0);
    arma::mat Arow(m*m,m*m*(p+1));
    Arow.fill(0);
    Arow.cols(m*m*p,m*m*(p+1)-1) = id2mat;
    for(int i = 0; i < p; i++)
    {
      Arow.cols(m*m*i,m*m*(i+1)-1) = -1*kron(idmat,phi.cols((p-i-1)*m,(p-i)*m-1));
    }
    for(int i = 0; i <= p; i++)
    {
      Amat.submat(i*m*m,i*m*m,(i+1)*m*m-1,(i+p+1)*m*m-1) = Arow;
    }
    arma::mat newA = Amat.submat(0,0,(p+1)*m*m-1,p*m*m-1);
    for(int i = 0; i <= p; i++)
    {
      for(int j = 0; j < p; j++)
      {
        newA.submat(i*m*m,j*m*m,(i+1)*m*m-1,(j+1)*m*m-1) = newA.submat(i*m*m,j*m*m,(i+1)*m*m-1,(j+1)*m*m-1) * Kmat;
      }
    }
    arma::mat Asubmat((p+1)*m*m,(p+1)*m*m);
    Asubmat.fill(0);
    Asubmat.cols(0,m*m-1) = Amat.cols(p*m*m,(p+1)*m*m-1);
    for(int i = 0; i < p; i++)
    {
      Asubmat.cols((i+1)*m*m,(i+2)*m*m-1) = Amat.cols((p+1+i)*m*m,(p+2+i)*m*m-1) + newA.cols((p-1-i)*m*m,(p-i)*m*m-1);
    }
    arma::mat Bmat(m*m*(q+1),m*m*(p+q+1));
    Bmat.fill(0);
    arma::mat Brow(m*m,m*m*(p+1));
    Brow.fill(0);
    Brow.cols(0,m*m-1) = id2mat;
    for(int i = 0; i < p; i++)
    {
      Brow.cols(m*m*(i+1),m*m*(i+2)-1) = -1*kron(phi.cols(i*m,(i+1)*m-1),idmat);
    }
    for(int i = 0; i <= q; i++)
    {
      Bmat.submat(i*m*m,i*m*m,(i+1)*m*m-1,(i+p+1)*m*m-1) = Brow;
    }
    arma::mat Bsubmat((q+1)*m*m,(q+1)*m*m);
    Bsubmat.fill(0);
    Bsubmat = Bmat.cols(0,(q+1)*m*m-1);
    arma::mat gam_mix = solve(Bsubmat,gamvec_ma);
    arma::mat gam_mixpad((p+1)*m*m,1);
    gam_mixpad.fill(0);
    if(p <= q)
    {
      gam_mixpad = gam_mix.rows(0,(p+1)*m*m-1);
    } else
    {
      gam_mixpad.rows(0,(q+1)*m*m-1) = gam_mix.rows(0,(q+1)*m*m-1);
    }
    arma::mat gam_arma = solve(Asubmat,gam_mixpad);
    gam_mix.reshape(m,m*(q+1));
    for(int i = 0; i <= q; i++)
    {
      gam_mixcube.slice(i) = gam_mix.cols(i*m,(i+1)*m-1);
    }
    gam_arma.reshape(m,m*(p+1));
    for(int i = 0; i <= p; i++)
    {
      gam_armacube.slice(i) = gam_arma.cols(i*m,(i+1)*m-1);
    }
    phi.reshape(m,m*p);    
  } else
  {
    gam_armacube.slice(0) = gam_ma.slice(0);
    gam_mixcube = gam_ma;
  }
    
  if(maxlag <= p)
  {
    gamma_final = gam_armacube.slices(0,maxlag);
  } else
  {
    gamma_final.slices(0,p) = gam_armacube;
    arma::cube gam_mixcubepad(m,m,maxlag+1);
    gam_mixcubepad.fill(0);
    if(maxlag > q)
    {
      gam_mixcubepad.slices(0,q) = gam_mixcube;
    } else
    {
      gam_mixcubepad = gam_mixcube.slices(0,maxlag);
    }
    for(int k = 0; k < (maxlag-p); k++)
    {
      int len = p+1+k;
      arma::mat acf = gam_mixcubepad.slice(len);
      if (p > 0)
      {
        arma::mat Cmat(m*p,m);
        Cmat.fill(0);
        for(int i = 0; i < p; i++)
        {
          Cmat.rows(i*m,(i+1)*m-1) = gamma_final.slice(len-i-1);
        }
        acf = acf + phi * Cmat;
      }
      gamma_final.slice(k+p+1) = acf;
    }
  }
 
  return gamma_final;
  
}
 