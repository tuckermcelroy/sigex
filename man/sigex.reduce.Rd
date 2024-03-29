% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.reduce.r
\name{sigex.reduce}
\alias{sigex.reduce}
\title{Determine a reduced rank model from a given fitted model}
\usage{
sigex.reduce(data.ts, param, mdl, thresh, modelflag)
}
\arguments{
\item{data.ts}{A T x N matrix ts object}

\item{param}{Model parameters entered into
a list object with an intuitive structure.}

\item{mdl}{The specified sigex model, a list object}

\item{thresh}{Lower bound on Schur complements}

\item{modelflag}{When TRUE, small Schur complements imply rank reduction
in the new model.  When modelFlag is FALSE, small Schur
complements are replaced by exp(thresh)}
}
\value{
list of mdl.red and par.red
   mdl.red: the new sigex model, a list object
	par.red: the new param for the new model
}
\description{
Determine a reduced rank model from a given fitted model
}
