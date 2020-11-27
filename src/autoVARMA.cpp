// Use the RcppArmadillo package
// Requires different header file from Rcpp.h
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cx_vec getEigenValues(arma::mat M)
{
  return arma::eig_gen(M);
}

// [[Rcpp::export]]
arma::cx_vec complexExp(arma::vec x)
{
  int len = x.size();
  arma::vec xreal(len);
  arma::vec ximag(len);
  for(int i = 0; i < len; i++)
  {
    xreal[i] = cos(x[i]);
    ximag[i] = sin(x[i]);
  }
  arma::cx_vec z = arma::cx_vec(xreal,ximag);
  return z;
}

// [[Rcpp::export]]
arma::cx_vec polymult(arma::cx_vec a, arma::cx_vec b)
{
  int lena = a.size();
  int lenb = b.size();
  int dimb = lenb -1;
  int lenc = lena + lenb -1;
  arma::cx_vec c(lenc);
  c.fill(0);
  for(int h = 0; h < lenc; h++)
  {
    int k = std::min(h,dimb);
    int l = std::max(h-lena+1,0);
    for(int j = l; j <= k; j++)
    {
      c[h] = c[h] + a[h-j]*b[j];
    }
  }

  return c;

}

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
List ar_adjoint(arma::cube poly_array)
{

  int p = poly_array.n_slices -1;
  int N = poly_array.n_cols;
  arma::mat poly_0 = poly_array.slice(0);
  arma::vec ar_poly(p*N+1);
  ar_poly.fill(0);
  ar_poly[0] = det(poly_0);
  int r = p*(N-1) + 1;
  arma::cube adj_array(N,N,r);
  arma::mat poly_inv = inv(poly_0);
  adj_array.slice(0) = ar_poly[0]*poly_inv;

  if(p > 0)
  {
    arma::mat poly_mat(N,N*p);
    for(int i = 0; i < p; i++)
    {
      poly_mat.cols(i*N,(i+1)*N-1) = poly_array.slice(i+1);
    }
    arma::mat poly_coefs = -1*poly_inv * poly_mat;
    arma::mat comp_mat;
    if(p == 1)
    {
      comp_mat = poly_coefs;
    } else
    {
      comp_mat = arma::eye<arma::mat>((p+1)*N,(p+1)*N);
      comp_mat.submat(0,N,N-1,N*(p+1)-1) = poly_coefs;
      comp_mat = comp_mat.rows(0,p*N-1);
      comp_mat = comp_mat.cols(N,(p+1)*N-1);
    }
    arma::cx_vec poly_evals = getEigenValues(comp_mat);
    arma::cx_vec ar_prod(1);
    ar_prod[0] = ar_poly[0];
    for(int j = 0; j < p*N; j++)
    {
      arma::cx_vec a_root(2);
      a_root[0] = 1.0;
      a_root[1] = -1.0*poly_evals[j];
      ar_prod = polymult(ar_prod,a_root);
    }
    ar_poly = arma::real(ar_prod);
    for(int j = 1; j < r; j++)
    {
      adj_array.slice(j) = ar_poly[j]*poly_inv;
      int l = std::min(p,j);
      for(int k = 0; k < l; k++)
      {
        adj_array.slice(j) = adj_array.slice(j) - poly_inv * poly_array.slice(k+1) * adj_array.slice(j-k-1);
      }
    }
  }

  return List::create(Named("adjoint") = adj_array,Named("det") = ar_poly);

}

// [[Rcpp::export]]
arma::cube auto_VARMA(arma::mat param,int p,int q,int ps,int qs,int season,int grid,int maxlag)
{

  //***********************************************************************
  //
  //	autoVARMA.cpp
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
  //	Purpose: computes autocovariances of SVARMA usng frequency domain
  //	Background: function computes autocovariances of SVARMA (p,q,ps,qs) from lag zero
  //		to maxlag, with array inputs phi and theta.  SVARMA equation:
  //	(1 - phi[1]B ... - phi[p]B^p) (1 - Phi[1]B^s ... - Phi[ps]B^{s*ps}) X_t =
  //    (1 + theta[1]B ... + theta[q]B^q) (1 + Theta[1]B^s ... + Theta[qs]B^{s*qs}) WN_t.
  //	Inputs:
  //    param: matrix of dimension m x (p+q+ps+qs+1),
  //        equals [ phi | theta | phiseas | thetaseas | sigma ]
  //		  phi: block matrix of dimension N x N*p of VAR coefficients
  //		  theta: block matrix of dimension N x N*q of VMA coefficients
  //      phiseas:  block matrix of dimension N x N*ps of SVAR coefficients
  //      thetaseas:  block matrix of dimension N x N*qs of SVMA coefficients
  //		  sigma: N x N covariance matrix of white noise
  //   grid: Riemann mesh size
  //  Notes: for pure VMA, leave out phi; for pure VAR, leave out theta.
  //	Outputs:
  //		autocovariances at lags 0 through maxlag, as array of dimension m x m x (maxlag+1)
  //
  //**************************************************************************

  int N = param.n_rows;
  arma::mat idmat = arma::eye<arma::mat>(N,N);
  arma::cube gamma(N,N,maxlag+1);
  gamma.fill(0);
  arma::cube phi_long(N,N,p+1);
  phi_long.slice(0) = idmat;
  if(p > 0)
  {
    for(int i = 0; i < p; i++)
    {
      phi_long.slice(i+1) = -1*param.cols(i*N,(i+1)*N-1);
    }
  }
  List out1 = ar_adjoint(phi_long);
  arma::cube phi_adjoint = out1["adjoint"];
  arma::vec phi_det = out1["det"];
  arma::cube theta_long(N,N,q+1);
  theta_long.slice(0) = idmat;
  if(q > 0)
  {
    for(int i = 0; i < q; i++)
    {
      theta_long.slice(i+1) = param.cols((p+i)*N,(p+i+1)*N-1);
    }
  }
  arma::cube phiseas_long(N,N,ps+1);
  phiseas_long.slice(0) = idmat;
  if(ps > 0)
  {
    for(int i = 0; i < ps; i++)
    {
      phiseas_long.slice(i+1) = -1*param.cols((p+q+i)*N,(p+q+i+1)*N-1);
    }
  }
  List out2 = ar_adjoint(phiseas_long);
  arma::cube phiseas_adjoint = out2["adjoint"];
  arma::vec phiseas_det = out2["det"];
  arma::cube thetaseas_long(N,N,qs+1);
  thetaseas_long.slice(0) = idmat;
  if(qs > 0)
  {
    for(int i = 0; i < qs; i++)
    {
      thetaseas_long.slice(i+1) = param.cols((p+q+ps+i)*N,(p+q+ps+i+1)*N-1);
    }
  }
  arma::mat sigma = param.cols((p+q+ps+qs)*N,(p+q+ps+qs+1)*N-1);
  int len = 2*grid+1;
  arma::vec freqs(len);
  freqs[grid] = 0;
  for(int k = 1; k <= grid; k++)
  {
    freqs[grid+k] = arma::datum::pi * k/grid;
    freqs[grid-k] = -1*arma::datum::pi * k/grid;
  }
  arma::cx_vec lambdas = complexExp(-1*freqs);

  for(int k = 0; k < len; k++)
  {
    arma::cx_mat phi_z(N,N);
    phi_z.fill(0);
    for(unsigned i = 0; i < phi_adjoint.n_slices; i++)
    {
      phi_z = phi_z + phi_adjoint.slice(i) * pow(lambdas,i)[k];
    }
    arma::cx_double phi_detz = 0;
    for(unsigned i = 0; i < phi_det.size(); i++)
    {
      phi_detz = phi_detz + phi_det[i] * pow(lambdas,i)[k];
    }
    arma::cx_mat theta_z(N,N);
    theta_z.fill(0);
    for(unsigned i = 0; i < theta_long.n_slices; i++)
    {
      theta_z = theta_z + theta_long.slice(i) * pow(lambdas,i)[k];
    }
    arma::cx_mat phiseas_z(N,N);
    phiseas_z.fill(0);
    for(unsigned i = 0; i < phiseas_adjoint.n_slices; i++)
    {
      phiseas_z = phiseas_z + phiseas_adjoint.slice(i) * pow(lambdas,i*season)[k];
    }
    arma::cx_double phiseas_detz = 0;
    for(unsigned i = 0; i < phiseas_det.size(); i++)
    {
      phiseas_detz = phiseas_detz + phiseas_det[i] * pow(lambdas,i*season)[k];
    }
    arma::cx_mat thetaseas_z(N,N);
    thetaseas_z.fill(0);
    for(unsigned i = 0; i < thetaseas_long.n_slices; i++)
    {
      thetaseas_z = thetaseas_z + thetaseas_long.slice(i) * pow(lambdas,i*season)[k];
    }

    double simpson = 0;
    double denom = 6*grid;
    if(k % 2 == 0)
    {
      simpson = 2/denom;
    }
    if(k % 2 == 1)
    {
      simpson = 4/denom;
    }
    if(k == 0)
    {
      simpson = 1/denom;
    }
    if(k == 2*grid)
    {
      simpson = 1/denom;
    }

    for(int h = 0; h <= maxlag; h++)
    {
      arma::cx_mat spec(N,N);
      spec = phiseas_z * phi_z * theta_z * thetaseas_z * sigma * thetaseas_z.t() * theta_z.t() * phi_z.t() * phiseas_z.t();
      spec = spec * pow(lambdas,-h)[k] * pow(abs(phi_detz),-2) * pow(abs(phiseas_detz),-2);
      gamma.slice(h) = gamma.slice(h) + arma::real(spec) * simpson;
    }

  }

  return gamma;

}
