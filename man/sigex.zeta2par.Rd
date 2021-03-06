% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigex.zeta2par.r
\name{sigex.zeta2par}
\alias{sigex.zeta2par}
\title{Transform zeta to param}
\usage{
sigex.zeta2par(zeta, mdlType, N)
}
\arguments{
\item{zeta}{This is the portion of the psi vector that
corresponds to t.s. models, such as cycles}

\item{mdlType}{This is a component of mdl (the specified sigex model),
cited as mdl[[2]]}

\item{N}{Cross-section dimension}
}
\value{
zeta.par: This is a portion of the full
		param list, corresponding to param[[3]]
}
\description{
Transform zeta to param
}
