% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.lpmse.r
\name{sigex.lpmse}
\alias{sigex.lpmse}
\title{Computes signal extraction mse arising from LP filtering
	of trend-cycle, with low pass cutoff, for trend and cycle each.}
\usage{
sigex.lpmse(param, mdl, trendcyclecomp, sigcomps, grid, cutoff)
}
\arguments{
\item{param}{model parameters entered into
a list object with an intuitive structure.}

\item{mdl}{The specified sigex model, a list object}

\item{trendcyclecomp}{The (single) index of the trend-cycle component}

\item{sigcomps}{Provides indices of a desired component that
is disjoint from trend-cycle, so that MSEs of
trend+sigcomps and cycle+sigcomps are computed.
 (Pass in sigcomps = NULL to just get trend and cycle MSEs.)}

\item{grid}{Desired number of frequencies for spectrum calculations}

\item{cutoff}{A number between 0 and pi, with all frequencies < cutoff preserved}
}
\value{
list object with mse.trend and mse.cycle
	mse.trend: N x N matrix, MSE of trend
	mse.cycle: N x N matrix, MSE of cycle
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
	this acgf and delta.  The error spectral density calculations
	are found in:
	"Casting Vector Time Series: Algorithms for Forecasting,
	Imputation, and Signal Extraction," McElroy (2018).
	param is the name for the model parameters entered into
	a list object with a more intuitive structure, whereas
	psi refers to a vector of real numbers containing all
	hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
}
