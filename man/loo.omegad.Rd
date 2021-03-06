% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loo.R
\name{loo.omegad}
\alias{loo.omegad}
\alias{loo}
\title{Compute leave-one-out (LOO) scores.}
\usage{
\method{loo}{omegad}(x, items = FALSE, ...)
}
\arguments{
\item{x}{omegad object.}

\item{items}{Logical (Default: FALSE). Whether to compute LOO per item, per observation (TRUE), or just per observation (FALSE).}

\item{...}{Arguments passed to \code{\link[loo]{loo}}.}
}
\value{
LOO object (if items = FALSE). List of LOO objects, one for each item (if items = TRUE).
}
\description{
LOO method for omegad objects.
}
\details{
LOO method for omegad objects.
}
\seealso{
\code{\link[loo]{loo}}
}
\author{
Stephen R. Martin
}
