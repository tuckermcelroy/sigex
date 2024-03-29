% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.extract.r
\name{sigex.extract}
\alias{sigex.extract}
\title{Computes signal extraction estimates with two standard errors}
\usage{
sigex.extract(data.ts, my.filter, mdl, param)
}
\arguments{
\item{data.ts}{A T x N matrix ts object}

\item{my.filter}{List object corresponding to the output of sigex.signal,
a list object of f.mat and v.mat.
f.mat: array of dimension c(T,N,T,N), where f.mat[,j,,k]
    is the signal extraction matrix that utilizes input series k
    to generate the signal estimate for series j.
v.mat: array of dimension c(T,N,T,N), where v.mat[,j,,k]
    is the error covariance matrix arising from input series k
    used to generate the signal estimate for series j.}

\item{mdl}{The specified sigex model, a list object}

\item{param}{The model parameters entered into a list object}
}
\value{
list object with extract, upp, and low
	extract: T x N matrix of the signal estimates
	upp: as extract, plus twice the standard error
	low: as extract, minus twice the standard error
}
\description{
Background:
	A sigex model consists of process x = sum y, for
	stochastic components y.  Each component process y_t
	is either stationary or is reduced to stationarity by
	application of a differencing polynomial delta(B), i.e.
		w_t = delta(B) y_t   is stationary.
	We have a model for each w_t process, and can compute its
	autocovariance function (acf), and denote its autocovariance
	generating function (acgf) via gamma_w (B).
	The signal extraction filter for y_t is determined from
	this acgf and delta.
	param is the name for the model parameters entered into
	a list object with a more intuitive structure, whereas
	psi refers to a vector of real numbers containing all
	hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
}
