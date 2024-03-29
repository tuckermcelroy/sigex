% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.getfrf.r
\name{sigex.getfrf}
\alias{sigex.getfrf}
\title{Computes frequency response function for desired signal and plots it}
\usage{
sigex.getfrf(data.ts, param, mdl, sigcomps, plotit = TRUE, grid)
}
\arguments{
\item{data.ts}{A T x N matrix ts object}

\item{param}{The model parameters entered into a list object}

\item{mdl}{The specified sigex model, a list object}

\item{sigcomps}{Indices of the latent components composing the signal}

\item{plotit}{Boolean flag for whether frf should be plotted;
only plots if N <= 3}

\item{grid}{Desired number of frequencies for spectrum calculations}
}
\value{
frf.comp:  array of dimension c(N,N,grid), with complex number entries
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
\details{
Notes: take grid >> len, else numerical issues arise
}
