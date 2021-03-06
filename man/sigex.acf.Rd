% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.acf.r
\name{sigex.acf}
\alias{sigex.acf}
\title{Compute the autocovariance function of a differenced latent component}
\usage{
sigex.acf(L.par, D.par, mdl, comp, mdlPar, delta, maxlag, freqdom = FALSE)
}
\arguments{
\item{L.par}{Unit lower triangular matrix in GCD of the component's
white noise covariance matrix.}

\item{D.par}{Vector of logged entries of diagonal matrix in GCD
of the component's white noise covariance matrix.}

\item{mdl}{The specified sigex model, a list object}

\item{comp}{Index of the latent component}

\item{mdlPar}{This is the portion of param
corresponding to mdl[[2]], cited as param[[3]]}

\item{delta}{Differencing polynomial (corresponds to eta(B) in Background)
written in format c(delta0,delta1,...,deltad)}

\item{maxlag}{Number of autocovariances required}

\item{freqdom}{A flag, indicating whether frequency domain acf routine should be used.}
}
\value{
x.acf: matrix of dimension N x N*maxlag, consisting of autocovariance
		matrices stacked horizontally
}
\description{
Compute the autocovariance function of a differenced latent component
}
